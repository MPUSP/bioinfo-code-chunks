import argparse
from Bio import Phylo  # pyright: ignore
from ete3 import Tree, TreeStyle, NCBITaxa  # pyright: ignore
from matplotlib import colormaps, colors  # pyright: ignore
from matplotlib import pyplot as plt  # pyright: ignore
from matplotlib.lines import Line2D  # pyright: ignore
from typing import Dict


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


def get_taxon(ncbi, node_name, rank):
    """
    Retrieve all taxonomic ranks for a given node name.
    :param ncbi: NCBI taxon database
    :param node_name: Node name in the phylogenetic tree
    :param rank: Target taxonomic rank
    """
    if not ncbi:
        return None
    if " | " in node_name:
        species = node_name.split(" | ")[1]
    else:
        species = node_name
    taxid = ncbi.get_name_translator([species]).get(species, [None])[0]
    if not taxid:
        return None
    lineage = ncbi.get_lineage(taxid)
    ranks = ncbi.get_rank(lineage)
    taxname = [
        ncbi.get_taxid_translator([tid])[tid]
        for tid, rnk in ranks.items()
        if rnk == rank
    ]
    return taxname[0] if taxname else None


# custom node function to draw different colors for nodes
def layout(node, color=None):
    if node.is_leaf() and ncbi:
        rank_value = get_taxon(ncbi, node.name, args.rank)
        color = rank_colors.get(str(rank_value), None)
        if color is not None:
            node.img_style["fgcolor"] = color
        else:
            node.img_style["fgcolor"] = "#CCCCCC"
        node.img_style["size"] = 5
    else:
        node.img_style["fgcolor"] = None
        node.img_style["size"] = 0
    node.img_style["vt_line_width"] = 2
    node.img_style["hz_line_width"] = 2


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate a phylogenetic tree visualization."
    )
    parser.add_argument("-i", "--input", required=True, help="Input phyloxml file")
    parser.add_argument("-o", "--output", required=True, help="Output SVG file")
    parser.add_argument("-s", "--style", default="", help="Tree style options")
    parser.add_argument(
        "-t",
        "--type",
        choices=["cladogram", "phylogram", "unrooted"],
        default="phylogram",
        help="Display type of tree, default: phylogram (correct scale of evol. distance)",
    )
    parser.add_argument(
        "-r",
        "--rank",
        default="genus",
        help="NCBI taxonomic rank to color nodes (default: genus)",
    )

    args = parser.parse_args()
    rank_colors: Dict[str, str] = {}

    # import phyloxml tree
    tree = Phylo.read(args.input, "phyloxml")

    # parse with ete3 lib
    ete_tree = build_tree_for_ete3(tree.root)

    # parse str input for styles to dict
    style_kwargs = {}
    if args.style:
        for item in args.style.split():
            if not item or "=" not in item:
                continue
            key, value = item.split("=", 1)
            style_kwargs[key] = value

    # get NCBI taxon database, downloaded only once
    try:
        ncbi = NCBITaxa()
    except Exception as e:
        print(f"Error initializing NCBITaxa: {e}")
        ncbi = None

    # get all taxons of choice (e.g. genus) for supplied species
    if ncbi:
        ranks = [get_taxon(ncbi, i.name, args.rank) for i in ete_tree.get_leaves()]
        unique_ranks = sorted({r for r in ranks if r})
        # generate colors of length ranks using cmap
        if unique_ranks:
            cmap = colormaps["Spectral"].resampled(len(unique_ranks))
            rank_colors.update(
                {rank: colors.to_hex(cmap(i)) for i, rank in enumerate(unique_ranks)}
            )
        # add a legend to the plot with key-color pairs from cmap
        legend_items = []
        for gen in unique_ranks:
            color = rank_colors.get(gen, None)
            if color is not None:
                legend_items.append(
                    Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="w",
                        markerfacecolor=color,
                        markersize=8,
                        label=gen,
                    )
                )
        if legend_items:
            legend_path = args.output.replace(".svg", "_legend.svg")
            legend_fig, legend_ax = plt.subplots(
                figsize=(4, max(2, len(legend_items) * 0.25))
            )
            legend_ax.legend(
                handles=legend_items,
                title="Genus",
                loc="center",
                frameon=True,
            )
            legend_ax.axis("off")
            legend_fig.tight_layout()
            legend_fig.savefig(legend_path, format="svg")
            plt.close(legend_fig)
            print(f"Legend written to: {legend_path}")

    # create tree style
    ts = TreeStyle()
    for k, v in style_kwargs.items():
        if v in ["True", "False"]:
            v = v == "True"
        elif v.isdigit():
            v = int(v)
        elif v.strip("-").replace(".", "", 1).isdigit():
            v = float(v)
        setattr(ts, k, v)
    ts.layout_fn = layout

    # tree type
    if args.type == "cladogram":
        # ts.force_topology = True  # another way to draw a cladogram
        ete_tree.dist = 0.1
        ete_tree.convert_to_ultrametric()
    elif args.type == "unrooted":
        ts.force_topology = False
        ete_tree.unroot()
        ete_tree.dist = 0
    elif args.type == "phylogram":
        ete_tree.dist = 0.1
    else:
        raise ValueError(f"Unknown tree type: {args.type}")

    # render tree
    tree_log = ete_tree.render(
        args.output,
        tree_style=ts,
        w=600,
        units="px",
        dpi=150,
    )
    print(f"Tree written to: {args.output}")
