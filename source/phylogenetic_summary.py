import argparse
from Bio import Phylo  # pyright: ignore
from ete3 import Tree, NCBITaxa  # pyright: ignore
from typing import Dict

parser = argparse.ArgumentParser(description="Summarize a phylogenetic tree")
parser.add_argument("-i", "--input", help="Input phyloxml file")
parser.add_argument("-o", "--output", help="Output SVG file")
args = parser.parse_args()


def build_tree_for_ete3(clade):
    """
    Convert a biopython tree to an ete3 tree
    :param clade: Biopython clade object
    """
    ete_tree = Tree()
    for child in clade.clades:
        child_tree = build_tree_for_ete3(child)
        ete_tree.add_child(child=child_tree, name=child.name, dist=child.branch_length)
    return ete_tree


def get_all_taxa(ncbi, node_name):
    """
    Retrieve all taxonomic ranks for a given node name.
    :param ncbi: NCBI taxon database
    :param node_name: Node name in the phylogenetic tree
    """
    if not ncbi:
        return {}
    if " | " in node_name:
        species = node_name.split(" | ")[1]
    else:
        species = node_name
    taxid = ncbi.get_name_translator([species]).get(species, [None])[0]
    if not taxid:
        return {}
    lineage = ncbi.get_lineage(taxid)
    ranks = ncbi.get_rank(lineage)
    taxa = {rnk: tid for tid, rnk in ranks.items()}
    return taxa if taxa else {}


def get_target_taxon(taxid, target_rank, ncbi):
    """
    For a given taxid and taxonomic rank (e.g. class 'bacilli'),
    retrieve the set of all neighboring taxa of the same rank
    :param taxid: Taxonomic ID
    :param target_rank: Target taxonomic rank
    :param ncbi: NCBI taxon database
    """
    lineage = ncbi.get_lineage(taxid)
    ranks = {rank: taxid for taxid, rank in ncbi.get_rank(lineage).items()}
    parent_rank = rank_order[rank_order.index(target_rank) - 1]
    try:
        parent_taxid = ranks[parent_rank]
    except KeyError:
        return []
    children = ncbi.get_descendant_taxa(parent_taxid)
    children_final = []
    for child in children:
        child_lineage = ncbi.get_lineage(child)
        ranks = ncbi.get_rank(child_lineage)
        children_final += [
            ncbi.get_taxid_translator([tid])[tid]
            for tid, rnk in ranks.items()
            if rnk == target_rank
        ]
    return list(set(children_final))


if __name__ == "__main__":

    # import phyloxml tree
    tree = Phylo.read(args.input, "phyloxml")

    # parse with ete3 lib
    ete_tree = build_tree_for_ete3(tree.root)

    # get NCBI taxon database, downloaded only once
    try:
        ncbi = NCBITaxa()
    except Exception as e:
        print(f"Error initializing NCBITaxa: {e}")
        raise e

    # loop through all leaf nodes and retrieve ranks as dict
    dict_taxa: Dict[str, Dict] = {}
    if ncbi:
        for leaf in ete_tree.get_leaves():
            dict_taxa[leaf.name] = get_all_taxa(ncbi, leaf.name)

    rank_order = [
        "domain",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]

    result = []
    for target_rank in rank_order[1:7]:
        # from the tree of known taxa, extract the target rank
        target_taxa = [
            v
            for item in dict_taxa.values()
            for k, v in item.items()
            if k == target_rank
        ]
        target_taxa = list(set(target_taxa))
        target_taxa_names = [
            ncbi.get_taxid_translator([tid])[tid] for tid in target_taxa
        ]
        # retrieve _all other_ taxa of the specified rank and taxid
        all_taxa = []
        for t in target_taxa:
            all_taxa += get_target_taxon(t, target_rank, ncbi)
        all_taxa = list(set(all_taxa))
        orphan_taxa = [t for t in target_taxa_names if t not in all_taxa]
        all_taxa = list(set(all_taxa) | set(orphan_taxa))
        neighboring_taxa = list(set(all_taxa) - set(target_taxa_names))
        result_formatted = "\t".join(
            [
                target_rank,
                str(len(target_taxa)),
                str(len(all_taxa)),
                str(len(neighboring_taxa)),
                "; ".join(target_taxa_names),
                "; ".join(all_taxa),
                "; ".join(neighboring_taxa),
            ]
        )
        print("Found for next rank:", result_formatted)
        result.append(result_formatted)

    # write summary to tsv file
    with open(args.output, "w") as out_tsv:
        out_tsv.write(
            "rank\ttarget\tall\tneighboring\ttarget_names\tall_names\tneighboring_names\n"
        )
        out_tsv.write("\n".join(result))
