#!/usr/bin/env python3
import argparse
import pandas as pd
from matplotlib import pyplot as plt


def annotate_stacked_bars(ax, fmt="{:.0f}"):
    """
    Annotate stacked bar chart with values on top of each segment.
    :param ax: selected axis
    :param fmt: format string for annotation
    """
    total_height = max(
        [
            rect.get_height()
            for container in ax.containers
            for rect in container
            if rect.get_height() > 0
        ]
    )
    for container in ax.containers:
        for rect in container:
            height = rect.get_height()
            if height <= 0:
                continue
            x = rect.get_x() + rect.get_width() / 2
            y = rect.get_y() + height + total_height / 30
            ax.text(
                x,
                y,
                fmt.format(height),
                ha="center",
                va="center",
                color="white",
                fontsize=8,
            )


if __name__ == "__main__":

    # parse arguments
    parser = argparse.ArgumentParser(
        description="Plot a summary of presence/absence data"
    )
    parser.add_argument("-i", "--input", required=True, help="Input summary TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output plot SVG file")
    args = parser.parse_args()

    # import tsv results files using pandas
    df = pd.read_csv(args.input, sep="\t")
    df["target_pc"] = df["target"] / df["all"] * 100
    df["neighbor_pc"] = df["neighboring"] / df["all"] * 100

    # plot a stacked barchart with pyplot using column 'rank' as names
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    df.loc[:, ["rank", "target", "neighboring"]].set_index("rank").plot(
        kind="bar", stacked=True, color=["red", "skyblue"], ax=axes[0]
    )
    df.loc[:, ["rank", "target_pc", "neighbor_pc"]].set_index("rank").plot(
        kind="bar", stacked=True, color=["red", "skyblue"], ax=axes[1]
    )
    axes[0].set_title("presence in all bacterial taxa", loc="left", fontsize=10)
    axes[1].set_title("presence in all bacterial taxa (%)", loc="left", fontsize=10)
    axes[0].legend(title="", fontsize=10, loc="upper right")
    axes[1].legend().remove()

    annotate_stacked_bars(axes[0])
    annotate_stacked_bars(axes[1], fmt="{:.1f}")
    plt.savefig(args.output)
    plt.close(fig)
    print(f"Plot saved to: {args.output}")
