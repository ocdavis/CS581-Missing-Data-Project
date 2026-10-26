import argparse
import os
from utils import load_config, sanitize_network
from pathlib import Path

CONFIG_HELP = r"""
The TOML configuration file must contain:

network = "/path/to/network/directories"
gene_trees = "/path/to/gene/tree/directories"
output = "/path/to/output/directory"
threads = 8
phylonet_path = "/path/to/PhyloNet.jar"
dry_run = false
gene_tree_method = "gene trees file name"
network_method = "network file name"
reticulations = 2
fixed = false

Example:

    network = "/data/networks"
    gene_trees = "/data/gene_trees"
    output = "/data/results"
    threads = 4
    phylonet_path = "/software/PhyloNet.jar"
    dry_run = false
    gene_tree_method = "iqtree_500-rooted.nwk"
    network_method = "camus.nwk"
    reticulations = 2
    fixed = false
"""

REQUIRED_KEYS = {
    "network",
    "gene_trees",
    "output",
    "threads",
    "phylonet_path",
    "dry_run",
    "gene_tree_method",
    "network_method",
    "reticulations",
    "fixed",
    "o_flag"
}


def create_nexus_file(path, gts, st, nret, threads, fixed, o_flag):
    g_name, s_name = (
        os.path.basename(gts).rsplit(".", 1)[0],
        os.path.basename(st).rsplit(".", 1)[0],
    )
    identifier = f'{g_name}_{s_name}_{nret}{"_fixed" if fixed else ""}'
    output_fname = os.path.join(path, f"phylonet_output_{identifier}.txt")
    nex_file = os.path.join(path, f"{identifier}.nex")
    with open(gts, "r") as f_gts, open(st, "r") as f_st, open(nex_file, "w") as f_nex:
        nex_contents = f"""#NEXUS

BEGIN NETWORKS;

Network st = {sanitize_network(f_st.read().strip())}

END;

BEGIN TREES;

"""
        for i, gt in enumerate(f_gts):
            nex_contents += f"Tree gt{i} = {gt.strip()}\n"
        nex_contents += f"""
END;

BEGIN PHYLONET;

InferNetwork_MPL (all) {nret} -s st {'-fs' if fixed else ''} {'-o' if o_flag else ''} -pl {threads} {output_fname};

END;
"""
        f_nex.write(nex_contents)
    return nex_file


def infer_all_phylonet(network_path, gene_tree_path, output_path, threads, phylonet_path, dry_run, gene_tree_method, network_method, reticulations, fixed, o_flag):
    stats_path = os.path.join(output_path, "stats.csv")
    for root, dirs, files in os.walk(network_path):
        sub_path = os.path.relpath(root, network_path)
        replicate_gene_tree_path = os.path.join(gene_tree_path, sub_path, gene_tree_method)
        if not os.path.exists(replicate_gene_tree_path):
            print(f"{replicate_gene_tree_path} does not exist, continuing.")
            continue
        replicate_output_path = os.path.join(output_path, sub_path)
        os.makedirs(replicate_output_path, exist_ok=True)
        for file in files:
            if file == network_method:
                nexus_file = create_nexus_file(replicate_output_path, replicate_gene_tree_path, os.path.join(root, file), reticulations, threads, fixed, o_flag)
                if not dry_run:
                    f_name = os.path.splitext(file)[0]
                    log_path = os.path.join(replicate_output_path, f"phylonet-{f_name}.log")
                    command = f'/usr/bin/time -a -f "{sub_path},{f_name},%e,%M" -o "{stats_path}"  java -jar {phylonet_path} {nexus_file} > {log_path}'
                    os.system(command)
                    Path("/u/quadros3/research/phylogeny/project_01/experiment_01/1a/scripts/sbatch/logfile.txt").unlink(missing_ok=True)
                break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to create nexus file and run PhyloNet MPL's InferNetwork_MPL for all replicates in a directory",
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
    infer_all_phylonet(
        config["network"],
        config["gene_trees"],
        config["output"],
        config["threads"],
        config["phylonet_path"],
        config["dry_run"],
        config["gene_tree_method"],
        config["network_method"],
        config["reticulations"],
        config["fixed"],
        config["o_flag"]
    )
