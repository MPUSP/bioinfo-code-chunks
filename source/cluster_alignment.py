import argparse
import re
import matplotlib.pyplot as plt  # pyright: ignore
from pathlib import Path
from typing import List
from Bio import AlignIO  # pyright: ignore
from Bio.Align import MultipleSeqAlignment  # pyright: ignore
from Bio import Phylo  # pyright: ignore
from Bio import SeqIO  # pyright: ignore
from Bio.Phylo.TreeConstruction import (  # pyright: ignore
    DistanceCalculator,
    DistanceMatrix,
    DistanceTreeConstructor,
)


def cluster_alignment(dm, args):
    """
    Cluster sequences based on a distance matrix and similarity threshold.
    :param dm: Distance matrix representing pairwise distances between sequences
    :param args: Command-line arguments containing clustering parameters
    """
    # compute similarity of sequences to each other
    identity = []
    total = len(dm)
    if args.similarity_metric == "average":
        for i in range(len(dm)):
            others = [1 - dm[i, j] for j in range(total) if j != i]
            identity.append(sum(others) / len(others) if others else 1.0)
    elif args.similarity_metric == "min":
        for i in range(total):
            row = [j for j in dm[i] if j != 0]
            identity.append(1 - min(row) if row else 1.0)
    else:
        raise ValueError(f"Unknown similarity metric: {args.similarity_metric}")
    #
    # sort by identity
    order = sorted(range(total), key=lambda idx: identity[idx], reverse=True)
    unassigned = set(range(total))
    clusters: list[list[int]] = []
    #
    # group sequences with identity above threshold
    for idx in order:
        if idx not in unassigned:
            continue
        members = [idx]
        for other in sorted(unassigned - {idx}):
            if 1 - dm[idx, other] >= args.threshold:
                members.append(other)
        for member in members:
            unassigned.discard(member)
        clusters.append(members)
    return clusters


def representative(cluster: List[int], dm: DistanceMatrix) -> int:
    """
    Get the most representative member of a cluster.
    :param cluster: List of sequence indices in the cluster
    :type cluster: List[int]
    :return: Index of the most representative sequence
    :rtype: int
    """
    if len(cluster) == 1:
        return cluster[0]
    best = cluster[0]
    best_score = -1.0
    for candidate in cluster:
        score = sum(1 - dm[candidate, other] for other in cluster if other != candidate)
        if len(cluster) > 1:
            score /= len(cluster) - 1
        if score > best_score:
            best_score = score
            best = candidate
    return best


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Cluster an alignment and write representative sequences"
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Aligned fasta file (MUSCLE output)"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Filtered fasta with representatives"
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        choices=["identity", "blastp", "blosum45", "blosum62"],
        default="identity",
        help="Distance calculation method",
    )
    parser.add_argument(
        "-s",
        "--similarity_metric",
        type=str,
        choices=["average", "min"],
        default="average",
        help="Similarity metric to use for grouping redundant sequences",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.9,
        help="Identity threshold for grouping redundant sequences",
    )
    parser.add_argument(
        "--tree-output",
        default=True,
        help="Write the phylogenetic tree image to svg",
    )

    # parse all arguments
    # testargs: ["-i", "results/blastp/non_spy/results_aligned.fasta", "-o", "results/blastp/non_spy/results_filtered.fasta"]
    args = parser.parse_args()

    # read input
    print(f"Reading alignment: {args.input}")
    alignment = AlignIO.read(args.input, "fasta")
    records = list(alignment)
    print(f"Aligned sequences: {len(records)}")

    # compute the distance matrix
    calculator = DistanceCalculator(args.method)
    dm = calculator.get_distance(alignment)

    # cluster the matrix by similarity
    clusters = cluster_alignment(dm, args)

    # get most representative member of a cluster
    representatives: List = []
    for index, cluster in enumerate(clusters, start=1):
        rep_idx = representative(cluster, dm)
        records[rep_idx].description += f" cluster_size={len(cluster)}"
        representatives.append(records[rep_idx])
        print(
            f"Cluster {index}: {len(cluster)} seqs -> representative {records[rep_idx].id}"
        )

    # build a tree ('upgma' or neighbor joining 'nj')
    constructor = DistanceTreeConstructor()
    alignment_filt = MultipleSeqAlignment(representatives)
    for seqrec in alignment_filt:
        species = re.search(r"\[(.*?)\]", seqrec.description)
        spec_name = f"{seqrec.name} | {species.group(1)}" if species else seqrec.name
        seqrec.id = spec_name
    dm_filt = calculator.get_distance(alignment_filt)
    tree = constructor.nj(dm_filt)
    out_prefix = Path(args.output).with_suffix("")

    # write and visualize tree
    Phylo.write(tree, out_prefix.with_suffix(".xml"), "phyloxml")
    print(f"Exported phylogenetic tree as: {out_prefix.with_suffix('.xml')}")
    if args.tree_output:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=ax, do_show=False)
        fig.savefig(str(out_prefix.with_suffix(".svg")), format="svg")
        plt.close(fig)
        print(f"Exported phylogenetic tree as: {out_prefix.with_suffix('.svg')}")

    # write the filtered alignment with representatives only
    SeqIO.write(representatives, out_prefix.with_suffix(".fasta"), "fasta")
    print(
        f"Filtered alignment with {len(representatives)} records written to: {out_prefix.with_suffix('.fasta')}"
    )

    # save the WP protein IDs to a text file
    with open(str(out_prefix.with_suffix(".txt")), "w") as out_stream:
        for wp in representatives:
            out_stream.write(f"{wp.id.split(' | ')[0]}.1\n")
    print(
        f"File with {len(representatives)} protein IDs written to: {out_prefix.with_suffix('.txt')}"
    )
