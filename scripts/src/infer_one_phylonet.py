import argparse
import os
from utils import load_config
from infer_all_phylonet import create_nexus_file
from pathlib import Path

CONFIG_HELP = r"""
The TOML configuration file must contain:

network = "/path/to/network.nwk"
gene_trees = "/path/to/gene/trees.nwk"
output = "/path/to/output/directory"
threads = 8
phylonet_path = "/path/to/PhyloNet.jar"
dry_run = false
reticulations = 2
fixed = false
o_flag = false

Unlike infer_all_phylonet.py, network and gene_trees point directly at a
single network file and a single gene trees file (no directories walked).

Example:

    network = "/data/networks/true_net-rooted.nwk"
    gene_trees = "/data/gene_trees/iqtree_500-rooted.nwk"
    output = "/data/results"
    threads = 4
    phylonet_path = "/software/PhyloNet.jar"
    dry_run = false
    reticulations = 1
    fixed = false
    o_flag = false
"""

REQUIRED_KEYS = {
    "network",
    "gene_trees",
    "output",
    "threads",
    "phylonet_path",
    "dry_run",
    "reticulations",
    "fixed",
    "o_flag",
}


def infer_one_phylonet(network_path, gene_tree_path, output_path, threads, phylonet_path, dry_run, reticulations, fixed, o_flag):
    if not os.path.isfile(network_path):
        print(f"{network_path} does not exist, exiting.")
        return
    if not os.path.isfile(gene_tree_path):
        print(f"{gene_tree_path} does not exist, exiting.")
        return
    os.makedirs(output_path, exist_ok=True)
    stats_path = os.path.join(output_path, "stats.csv")
    nexus_file = create_nexus_file(output_path, gene_tree_path, network_path, reticulations, threads, fixed, o_flag)
    if not dry_run:
        f_name = os.path.splitext(os.path.basename(network_path))[0]
        log_path = os.path.join(output_path, f"phylonet-{f_name}.log")
        command = f'/usr/bin/time -a -f "{f_name},%e,%M" -o "{stats_path}"  java -jar {phylonet_path} {nexus_file} > {log_path}'
        os.system(command)
        Path("/u/quadros3/research/phylogeny/project_01/experiment_01/1a/scripts/sbatch/logfile.txt").unlink(missing_ok=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to create nexus file and run PhyloNet MPL's InferNetwork_MPL for a single network",
        epilog=CONFIG_HELP,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-c",
        "--config",
        required=True,
        help="path to the TOML configuration file",
    )

    args = parser.parse_args()
    config = load_config(args.config, REQUIRED_KEYS)

    print(config)
    infer_one_phylonet(
        config["network"],
        config["gene_trees"],
        config["output"],
        config["threads"],
        config["phylonet_path"],
        config["dry_run"],
        config["reticulations"],
        config["fixed"],
        config["o_flag"],
    )
