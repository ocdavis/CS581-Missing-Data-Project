"""Script to create Nexus file and run PhyloNet MPL's Cmpnets
Fall 2026, Sathvik Quadros
"""

import argparse
import os
from run_cmpnets import create_nexus_file
from utils import load_config

CONFIG_HELP = r"""
The TOML configuration file must contain:

true_network = "/path/to/true/network/directories"
computed_network = "/path/to/computed/network/directories"
output = "/path/to/output/directory"
mode = "characterization mode"
phylonet_path = "/path/to/PhyloNet.jar"
dry_run = false
true_net_name = "true_net.nwk"
network_methods = ["list", "of", "network", "filenames"]

Example:

    true_network = "/data/true_networks"
    computed_network = "/data/computed_networks"
    output = "/data/results"
    mode = "cluster"
    phylonet_path = "/software/PhyloNet.jar"
    dry_run = false
    true_net_name = "true_net.nwk"
    network_methods = ["camus", "phylonet-", "phylonet_ft-", "reroot-snaq"]
"""

REQUIRED_KEYS = {
    "true_network",
    "computed_network",
    "output",
    "mode",
    "phylonet_path",
    "dry_run",
    "true_net_name",
    "network_methods"
}


def run_all_cmpnets(tnet_path, cnet_path, output_path, mode, phylonet_path, dry_run, true_net_name, network_methods):
    stats_path = os.path.join(output_path, "stats.csv")
    with open(stats_path, "w") as f:    
        f.write("network_file,output_file,replicate,method_name,time,memory")
    for root, dirs, files in os.walk(cnet_path):
        sub_path = os.path.relpath(root, cnet_path)
        replicate_tnet_path = os.path.join(tnet_path, sub_path, true_net_name)
        if not os.path.exists(replicate_tnet_path):
            print(f"{replicate_tnet_path} does not exist, continuing.")
            continue
        replicate_output_path = os.path.join(output_path, sub_path)
        os.makedirs(replicate_output_path, exist_ok=True)
        for file in files:
            if any(method in file for method in network_methods) and ".nwk" in file:
                replicate_cnet_path = os.path.join(root, file)
                # maybe change how the output file is being obtained for the purposes of the stats.csv
                (nexus_file, replicate_output_file_path) = create_nexus_file(replicate_output_path, replicate_tnet_path, replicate_cnet_path, mode) 
                if not dry_run:
                    f_name = os.path.splitext(file)[0]
                    cmpnets_log_path = os.path.join(replicate_output_path, f"cmpnet_{f_name}.log")
                    print(cmpnets_log_path)
                    command = f'/usr/bin/time -a -f "{replicate_cnet_path},{replicate_output_file_path},{sub_path},{f_name},%e,%M" -o "{stats_path}"  java -jar {phylonet_path} {nexus_file} > {cmpnets_log_path}'
                    os.system(command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to create nexus file and run PhyloNet MPL's Cmpnets for all replicates in a directory",
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
    run_all_cmpnets(
        config["true_network"],
        config["computed_network"],
        config["output"],
        config["mode"],
        config["phylonet_path"],
        config["dry_run"],
        config["true_net_name"],
        config["network_methods"]
    )
