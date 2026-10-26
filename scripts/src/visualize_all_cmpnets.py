"""Script to visualize all comparison networks
Fall 2026, Sathvik Quadros
"""
import argparse
import os
import pandas as pd
from utils import load_config, cmpnets_accumulate_records
import seaborn
import matplotlib.pyplot as plt


CONFIG_HELP = r"""
The TOML configuration file must contain:

input_path = "/path/to/input/file"
output_path = "/path/to/output/directory"

Example:

    input_path = "/u/quadros3/research/phylogeny/project_01/experiment_01/1a/output_cmpnets/n15/cmpnets_dataframe.csv"
    output_path = "/u/quadros3/research/phylogeny/project_01/experiment_01/1a/output_cmpnets"
"""

REQUIRED_KEYS = {
    "input_path",
    "output_path",
}


def visualize_all_cmpnets(input_path, output_path):
    """Generate boxplots and CSV from comparison network results."""
    df = pd.read_csv(input_path, dtype={"replicate": str})

    # Figure 1: False Negative
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    seaborn.boxplot(
        data=df,
        x="method_label",
        y="false_negative",
        ax=ax1
    )
    ax1.set_title("False Negative Distribution")
    ax1.set_xlabel("Method")
    ax1.set_ylabel("False Negative Rate")
    ax1.tick_params(axis="x", rotation=20)
    fig1.tight_layout()

    # Figure 2: False Positive
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    seaborn.boxplot(
        data=df,
        x="method_label",
        y="false_positive",
        ax=ax2
    )
    ax2.set_title("False Positive Distribution")
    ax2.set_xlabel("Method")
    ax2.set_ylabel("False Positive Rate")
    ax2.tick_params(axis="x", rotation=20)
    fig2.tight_layout()

    # Figure 3: False Average
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    seaborn.boxplot(
        data=df,
        x="method_label",
        y="false_average",
        ax=ax3
    )
    ax3.set_title("Average of FP & FN Distribution")
    ax3.set_xlabel("Method")
    ax3.set_ylabel("Average Error Rate")
    ax3.tick_params(axis="x", rotation=20)
    fig3.tight_layout()

    # Save outputs
    os.makedirs(output_path, exist_ok=True)
    fig1.savefig(os.path.join(output_path, "cmpnets_fn_boxplot.png"))
    fig2.savefig(os.path.join(output_path, "cmpnets_fp_boxplot.png"))
    fig3.savefig(os.path.join(output_path, "cmpnets_fa_boxplot.png"))

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

    visualize_all_cmpnets(
        config["input_path"],
        config["output_path"],
    )
