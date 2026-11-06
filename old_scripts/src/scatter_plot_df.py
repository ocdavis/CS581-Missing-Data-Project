"""Script to create Nexus file and run PhyloNet MPL's CalGTProb
Fall 2026, Sathvik Quadros
"""
import pandas as pd
import argparse
import os

def scatter_plot(calgtprob_input, cmpnets_input, output_path):
    df_cal = pd.read_csv(calgtprob_input, dtype={"replicate": str})
    df_cmp = pd.read_csv(cmpnets_input, dtype={"replicate": str})

    df_cmb = df_cal.merge(df_cmp, on=['replicate', 'method_label'])

    y_axis = 'false_positive'
    print(df_cmb)
    plot = df_cmb.plot.scatter(x='mpl_score', y=y_axis)
    fig = plot.get_figure()
    fig.savefig(os.path.join(output_path, f"output_{y_axis}.png"))
    # df_cmb.to_csv(os.path.join(output_path, f"comb.csv")) #, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to create table and rank cmpnet scores for different methods",
    )
    parser.add_argument(
        "-i1", "--input-path1", type=str, required=True, help="path to ranked calgtprob"
    )
    parser.add_argument(
        "-i2", "--input-path2", type=str, required=True, help="path to ranked cmpnets"
    )
    parser.add_argument(
        "-o", "--output-path", type=str, required=True, help="output directory"
    )
    args = parser.parse_args()
    scatter_plot(
        args.input_path1,
        args.input_path2,
        args.output_path,
    )
