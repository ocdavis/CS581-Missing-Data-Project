"""Script to create Nexus file and run PhyloNet MPL's Cmpnets
Modified from https://doi.org/10.13012/B2IDB-6892704_V1} by James Willson
Fall 2026, Sathvik Quadros
"""

import argparse
import os
from src.utils import remove_branch_lengths


def create_nexus_file(path, tnet, cnet, mode, file_name=None):
    if file_name is None:
        file_name = os.path.splitext(os.path.basename(cnet))[0]
    output_fname = os.path.join(path, f"cmpnets_output_{file_name}.txt")
    nex_file = os.path.join(path, f"cmpnets_{file_name}.nex")
    with open(tnet, "r") as f_tnet, open(cnet, "r") as f_cnet, open(nex_file, "w") as f_nex:
        nex_contents = f"""#NEXUS

BEGIN NETWORKS;

Network tnet = {remove_branch_lengths(f_tnet.read().strip())}
Network cnet = {remove_branch_lengths(f_cnet.read().strip())}

END;

BEGIN PHYLONET;

Cmpnets tnet cnet -m {mode} {output_fname};

END;
"""
        f_nex.write(nex_contents)
    return (nex_file, output_fname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="script to create nexus file and run PhyloNet MPL's Cmpnets",
    )
    parser.add_argument(
        "-tn", "--true-network", type=str, required=True, help="true network file"
    )
    parser.add_argument(
        "-cn", "--computed-network", type=str, required=True, help="computed network file"
    )
    parser.add_argument(
        "-m", "--mode", type=str, required=True, help="Characterization mode, read phylonet documentation: https://phylogenomics.rice.edu/html/commands/Cmpnets.html"
    )
    parser.add_argument(
        "-p", "--phylonet-path", type=str, required=True, help="path to PhyloNet.jar"
    )
    parser.add_argument("-l", "--log", type=str, help="path to save log file")
    parser.add_argument(
        "-d",
        "--dry-run",
        action="store_true",
        help="dry run -- don't run phylonet",
    )
    args = parser.parse_args()
    path = os.path.dirname(args.true_network)
    nexus_file = create_nexus_file(
        path,
        args.true_network,
        args.computed_network,
        args.mode,
    )
    if not args.dry_run:
        command = f"java -jar {args.phylonet_path} {nexus_file}"
        if args.log is not None:
            command += f" > {args.log}"
        os.system(command)
