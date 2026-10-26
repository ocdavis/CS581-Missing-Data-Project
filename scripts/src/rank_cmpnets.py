"""Script to create Nexus file and run PhyloNet MPL's CalGTProb
Fall 2026, Sathvik Quadros
"""
import pandas as pd
import argparse
import os

def rank_cmpnets(input_path, output_path):
    df = pd.read_csv(input_path, dtype={"replicate": str})
    categories = ["false_positive", "false_negative", "false_average"]
    for category in categories:
        result = (
            df.assign(
                placement=df.groupby("replicate")[category]
                .rank(ascending=True)
            )
            .pivot(index="replicate", columns="method_label", values="placement")
            .reset_index()
        )

        result = result[result.notna().all(axis=1)]


        methods = df["method_label"].unique()
        mean = result[methods].mean()
        mean.name = "mean"

        median = result[methods].median()
        median.name = "median"

        result = pd.concat([result, mean.to_frame().T])
        result = pd.concat([result, median.to_frame().T])

        print(result)
        result.to_csv(os.path.join(output_path, f"cmpnets_{category}.csv")) #, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to create table and rank cmpnet scores for different methods",
    )
    parser.add_argument(
        "-i", "--input-path", type=str, required=True, help="path to dataframe csv containing cmpnets scores"
    )
    parser.add_argument(
        "-o", "--output-path", type=str, required=True, help="output directory"
    )
    args = parser.parse_args()
    rank_cmpnets(
        args.input_path,
        args.output_path,
    )
