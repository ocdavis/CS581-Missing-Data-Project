import argparse
import os
import pandas as pd

def extract_camus_nwk(replicates_path, reticulations=1):
    for root, dirs, files in os.walk(replicates_path):
        for file in files:
            if "camus" in file and ".csv" in file:
                f_name = os.path.splitext(file)[0]
                nwk_path = os.path.join(root, f"{f_name}.nwk")
                df = pd.read_csv(os.path.join(root,file))
                newick = df.loc[df['Number of Branches']  == reticulations, "Extended Newick"].values[0]
                with open(nwk_path, "w") as f_nwk:
                    f_nwk.write(newick)
                break

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to create nexus file and run PhyloNet MPL's CalGTProb for all replicates in a directory",
    )
    parser.add_argument(
        "-n", "--network", type=str, required=True, help="path to parent of replicate directories containing computed networks"
    )
    args = parser.parse_args()
    extract_camus_nwk(
        args.network,
    )