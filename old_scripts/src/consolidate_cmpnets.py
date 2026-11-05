"""Script to visualize all comparison networks
Fall 2026, Sathvik Quadros
"""
import argparse
import os
import pandas as pd
from utils import load_config, cmpnets_consolidate_records, merge_records_with_stats


CONFIG_HELP = r"""
The TOML configuration file must contain:

input_path = "/path/to/input/directories"
output_path = "/path/to/output/directory"
file_names = ["filename1.txt", "filename2.txt", ...]

[method_labels]
filename1.txt = "Method Label 1"
filename2.txt = "Method Label 2"

Example:

    input_path = "/u/quadros3/research/phylogeny/project_01/experiment_01/1a/output_cmpnets"
    output_path = "/u/quadros3/research/phylogeny/project_01/experiment_01/1a/output_cmpnets"
    
    file_names = [
        "cmpnets_output_camus_f05_astral_iqtree.txt",
        "cmpnets_output_phylonet_ft-astral-iqtree-1.txt",
        "cmpnets_output_phylonet-astral-iqtree-1.txt",
        "cmpnets_output_reroot-snaq-astral-iqtree-1.txt"
    ]
    
    [method_labels]
    cmpnets_output_camus_f05_astral_iqtree.txt = "CAMUS"
    cmpnets_output_phylonet_ft-astral-iqtree-1.txt = "PhyloNET-FT"
    cmpnets_output_phylonet-astral-iqtree-1.txt = "PhyloNET"
    cmpnets_output_reroot-snaq-astral-iqtree-1.txt = "SNaQ"
"""

REQUIRED_KEYS = {
    "input_path",
    "output_path",
    "file_names",
    "method_labels"
}


def consolidate_cmpnets(input_path, output_path, file_names, method_labels):
    """Generate boxplots and CSV from comparison network results."""

    df_records = cmpnets_consolidate_records(input_path, file_names, method_labels)
    
    df_stats = pd.read_csv(os.path.join(input_path, "stats.csv"), dtype={"replicate": str})

    df_merged = merge_records_with_stats(df_records, df_stats)

    # Save outputs
    os.makedirs(output_path, exist_ok=True)
    df_merged.to_csv(os.path.join(output_path, "cmpnets_dataframe.csv"), index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to visualize comparing all networks",
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

    consolidate_cmpnets(
        config["input_path"],
        config["output_path"],
        config["file_names"],
        config["method_labels"]
    )
