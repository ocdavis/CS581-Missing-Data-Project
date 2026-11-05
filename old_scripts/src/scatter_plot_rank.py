"""Script to create Nexus file and run PhyloNet MPL's CalGTProb
Fall 2026, Sathvik Quadros
"""
import pandas as pd
import tomllib
import os
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from utils import load_config

REQUIRED_KEYS = {"calgtprob_input", "cmpnets_input", "output_path", "methods"}

def scatter_plot(calgtprob_input, cmpnets_input, output_path, methods):
    # Read the files
    m_df = pd.read_csv(calgtprob_input, dtype={"replicate": str})
    f_df = pd.read_csv(cmpnets_input, dtype={"replicate": str})

    m_df = m_df[m_df.notna().all(axis=1)]
    f_df = f_df[f_df.notna().all(axis=1)]

    # Ensure replicate IDs have the same format in both files
    m_df["replicate"] = m_df["replicate"].astype(str).str.zfill(2)
    f_df["replicate"] = f_df["replicate"].astype(str).str.zfill(2)


    # Convert both data frames from wide format to long format
    m_long = m_df.melt(
        id_vars="replicate",
        value_vars=methods,
        var_name="method",
        value_name="mpl_score"
    )

    f_long = f_df.melt(
        id_vars="replicate",
        value_vars=methods,
        var_name="method",
        value_name="FA_score"
    )

    # Match each M score with its corresponding F score
    combined = pd.merge(
        m_long,
        f_long,
        on=["replicate", "method"],
        how="inner"
    )

    print(combined)

    fig1, ax1 = plt.subplots(figsize=(10, 6))

    sns.stripplot(
        data=combined,
        x="mpl_score",
        y="FA_score",
        jitter=0.1,  # Add perturbation in order to make data more visible
        ax=ax1
    )

    ax1.set_title("MPL vs False Average ranked")
    ax1.set_xlabel("MPL Score Rank")
    ax1.set_ylabel("Average Error Rate Rank")

    fig1.tight_layout()

    os.makedirs(output_path, exist_ok=True)

    fig1.savefig(os.path.join(output_path, "mpl_vs_f.png"))
    combined.to_csv(os.path.join(output_path, "combined_mpl_fa.csv"), index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to visualize comparing all networks",
    )

    parser.add_argument(
        "-c",
        "--config",
        required=True,
        help="path to the TOML configuration file",
    )

    args = parser.parse_args()
    config = load_config(args.config, REQUIRED_KEYS)

    scatter_plot(
        config["calgtprob_input"],
        config["cmpnets_input"],
        config["output_path"],
        config["methods"]
    )
