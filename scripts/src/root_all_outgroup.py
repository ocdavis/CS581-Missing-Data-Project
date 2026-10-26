"""Script to root all unrooted gene trees in a directory and subdirs
Modified from https://doi.org/10.13012/B2IDB-6892704_V1} by James Willson
Fall 2026, Sathvik Quadros
"""
import treeswift as ts
import warnings
import argparse
import os

COMPUTED_GT_FNAMES = ["iqtree_500.nwk", "g_500.nwk", "g_true.nwk"]

def reroot_on_edge(tree, node):
    # suppress reroot is maybe buggy warning from treeswift
    with warnings.catch_warnings(action="ignore"):
        if not node.is_root():
            if (
                not hasattr(node, "edge_length")
                or node.edge_length is None
                or node.edge_length == 0
            ):
                node.edge_length = 1
            tree.reroot(node, length=node.edge_length / 2)

def root_outgroup(tree):
    if "OUT" in [c.label for c in tree.root.children]:
        return tree
    outgroup = [l for l in tree.traverse_leaves() if l.label == "OUT"][0]
    outgroup.edge_length = 1
    reroot_on_edge(tree, outgroup)
    tree.is_rooted = False
    tree.suppress_unifurcations()
    return tree

def root_single_outgroup(in_file, out_file):
    with open(in_file, "r") as fi, open(out_file, "w") as fo:
        for line in fi:
            fo.write(f"{root_outgroup(ts.read_tree_newick(line)).newick()}\n")

def root_all_outgroup(input_path, output_path):
    for root, dirs, files in os.walk(input_path):
        sub_path = os.path.relpath(root, input_path)
        for file in files:
            if any(method == file for method in COMPUTED_GT_FNAMES):
                file_read_path = os.path.join(root, file)
                parts = file.split('.')
                out_file_name = f"{parts[0]}-rooted.{parts[1]}"
                file_write_dir = os.path.join(output_path,sub_path)
                os.makedirs(file_write_dir,exist_ok=True)
                file_write_path = os.path.join(file_write_dir,out_file_name)
                root_single_outgroup(file_read_path, file_write_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to visualize all calgtprob mpl scores",
    )
    parser.add_argument(
        "-i", "--input-path", type=str, required=True, help="path to directory containing unrooted gene trees"
    )
    parser.add_argument(
        "-o", "--output-path", type=str, required=True, help="output directory path"
    )
    args = parser.parse_args()
    root_all_outgroup(
        args.input_path,
        args.output_path,
    )
