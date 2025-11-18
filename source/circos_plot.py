#!/usr/bin/env python3
"""
Generate Circos plots from genome FASTA and GFF annotation files.
Requires: pycircos, biopython, matplotlib
Install: pip install python-circos biopython matplotlib
"""

import sys
from Bio import SeqIO
import pycircos
import numpy as np
from collections import defaultdict


def parse_fasta(fasta_file):
    """Parse FASTA file and return sequence records."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = record
    return sequences


def parse_gff(gff_file, feature_types=None):
    """parse GFF file and extract features by type."""
    if feature_types is None:
        feature_types = ["CDS", "gene", "rRNA", "tRNA"]
    features = defaultdict(list)
    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attributes = parts
            if ftype in feature_types:
                features[seqid].append(
                    {
                        "type": ftype,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand,
                    }
                )
    return features


def create_circos_plot(fasta_file, gff_file, output_file="circos_plot.png"):
    """Create Circos plot with genome and annotations."""

    # parse input files
    print(f"parsing FASTA: {fasta_file}")
    sequences = parse_fasta(fasta_file)

    print(f"parsing GFF: {gff_file}")
    features = parse_gff(gff_file)

    # check that all seq_ids from annotation are present in seq
    seq_ids_fasta = list(sequences.keys())
    seq_ids_gff = list(features.keys())
    if not all(i in seq_ids_fasta for i in seq_ids_gff):
        raise ValueError("Some sequence IDs in GFF not found in FASTA.")

    # initialize Circos
    circle = pycircos.Gcircle(figsize=(10, 10))

    # add arcs for each sequence/chromosome
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
    arc_data = {}

    for idx, (seqid, record) in enumerate(sequences.items()):
        arc = pycircos.Garc(
            arc_id=seqid,
            size=len(record.seq),
            interspace=2,
            raxis_range=(900, 950),
            facecolor=colors[idx % len(colors)],
            edgecolor="black",
            linewidth=0.5,
            label=seqid,
            record=record,
        )
        circle.add_garc(arc)
        arc_data[seqid] = arc
        print(f"added arc: {seqid} (length: {len(record.seq):,} bp)")

    circle.set_garcs()

    # add tick labels for each arc
    for seqid in circle.garc_dict:
        ticks = range(0, round(arc_data[seqid].size / 1000), 100)
        circle.tickplot(
            seqid,
            raxis_range=(960, 985),
            tickinterval=100000,
            ticklabels=ticks,
            ticklabeldirection="outer",
        )

    # add features from GFF
    feature_colors = {
        "CDS": "#377eb8",
        "gene": "#4daf4a",
        "rRNA": "#e41a1c",
        "tRNA": "#ff7f00",
    }

    for seqid in features:
        if seqid not in arc_data:
            continue

        for feat in features[seqid]:
            circle.barplot(
                seqid,
                data=[1],
                positions=[feat["start"] - 1],
                width=[feat["end"] - feat["start"] + 1],
                raxis_range=[850, 890] if feat["strand"] == "+" else [810, 850],
                facecolor=[feature_colors.get(feat["type"], "#999999")],
            )

    # calculate GC content and GC skew
    print("calculating GC content and GC skew...")
    window = 1000
    for seqid, record in sequences.items():
        if seqid not in arc_data:
            continue
        gc_values = []

        # calculate GC skew for the entire sequence
        gc_skew = pycircos.Garc(
            arc_id=record.id, record=record, size=len(record.seq)
        ).calc_nnskew(n1="G", n2="C", window_size=window)
        gc_skew_positions = np.array(
            [(i * window + window // 2) for i in range(len(gc_skew))]
        )
        gc_skew_plus = np.clip(gc_skew, 0, None)
        gc_skew_minus = np.clip(gc_skew, None, 0)

        for i in range(0, len(record.seq), window):
            subseq = record.seq[i : i + window]
            if len(subseq) > 0:
                gc = (
                    subseq.count("G")
                    + subseq.count("C")
                    + subseq.count("g")
                    + subseq.count("c")
                ) / len(subseq)
                gc_values.append([i + window // 2, i + window, gc])

        if gc_values:
            min_skew = min(gc_skew_minus)
            max_skew = max(gc_skew_plus) 
            if min_skew == max_skew == 0:
                skew_rlim = (-0.1, 0.1)
            else:
                skew_rlim = (min_skew, max_skew)
            circle.lineplot(
                record.id,
                data=[i[2] for i in gc_values],
                positions=[i[0] for i in gc_values],
                raxis_range=(700, 800),
                rlim=(0.3, 0.7),
                linecolor="#984ea3",
                linewidth=1,
            )
            circle.fillplot(
                record.id,
                data=gc_skew_plus,
                positions=gc_skew_positions,
                raxis_range=(600, 700),
                base_value=0,
                rlim=skew_rlim,
                facecolor="#984ea3",
            )
            circle.fillplot(
                record.id,
                data=gc_skew_minus,
                positions=gc_skew_positions,
                raxis_range=(600, 700),
                base_value=0,
                rlim=skew_rlim,
                facecolor="#642f6c",
            )

    # save figure
    print(f"generating plot: {output_file}")
    circle.figure.savefig(output_file, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "usage: python circos_plot.py <genome.fasta> <annotation.gff> [output.png]"
        )
        sys.exit(1)

    fasta_file = sys.argv[1]
    gff_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else "circos_plot.png"
    create_circos_plot(fasta_file, gff_file, output_file)
