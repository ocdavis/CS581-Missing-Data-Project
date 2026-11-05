"""Script to compare same-rank inferred networks using Cmpnets
Fall 2026, Sathvik Quadros
"""

import argparse
import itertools
import math
import os
import re

import pandas as pd

from utils import load_config
from run_cmpnets import create_nexus_file
from utils import get_fn_fp_fa


CONFIG_HELP = r"""
The TOML configuration file must contain:

input_csv = "/path/to/cmpnets_dataframe.csv"
ranks_csv = "/path/to/cmpnets_false_average.csv"
output_path = "/path/to/output/directory"
phylonet_path = "/path/to/PhyloNet.jar"
mode = "cluster"
dry_run = false

Example:

    input_csv = "/u/quadros3/research/phylogeny/project_01/experiment_01/1a/cmpnets_dataframe.csv"
    ranks_csv = "/u/quadros3/research/phylogeny/project_01/experiment_01/1a/cmpnets_false_average.csv"
    output_path = "/u/quadros3/research/phylogeny/project_01/experiment_01/1a/same_rank_checks"
    phylonet_path = "/u/quadros3/software/PhyloNet.jar"
    mode = "cluster"
    dry_run = false
"""


REQUIRED_KEYS = {
    "input_csv",
    "ranks_csv",
    "output_path",
    "phylonet_path",
    "mode",
    "dry_run",
}


def slugify(value):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_") or "method"


def is_close(a, b):
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=1e-12)


def find(parent, method):
    while parent[method] != method:
        parent[method] = parent[parent[method]]
        method = parent[method]
    return method


def union(parent, method_a, method_b):
    root_a = find(parent, method_a)
    root_b = find(parent, method_b)

    if root_a == root_b:
        return root_a

    # Keep the lexicographically smaller name as the representative.
    if root_a < root_b:
        parent[root_b] = root_a
        return root_a

    parent[root_a] = root_b
    return root_b


def compare_same_rank_networks(input_csv, ranks_csv, output_path, phylonet_path, mode, dry_run):

    df = pd.read_csv(input_csv, dtype={"replicate": str})
    df["replicate"] = df["replicate"].astype(str).str.strip()
    df["method_label"] = df["method_label"].astype(str).str.strip()
    df = df.drop_duplicates(subset=["replicate", "method_label"], keep="first")

    row_lookup = df.set_index(["replicate", "method_label"]).to_dict(orient="index")

    ranks = pd.read_csv(ranks_csv, index_col=0, dtype={"replicate": str})
    ranks = ranks[~ranks.index.astype(str).str.strip().str.lower().isin(["mean", "median"])]
    ranks = ranks.dropna(subset=["replicate"])
    ranks["replicate"] = ranks["replicate"].astype(str).str.strip()
    ranks.columns = [str(c).strip() for c in ranks.columns]

    methods = [c for c in ranks.columns if c != "replicate"]

    os.makedirs(output_path, exist_ok=True)

    comparison_records = []
    group_records = []

    for _, rank_row in ranks.iterrows():
        replicate = rank_row["replicate"]

        parent = {method: method for method in methods}
        checked_representative_pairs = {}

        rank_groups = {}
        for method in methods:
            rank_value = rank_row[method]
            if pd.notna(rank_value):
                rank_key = round(float(rank_value), 12)
                rank_groups.setdefault(rank_key, []).append(method)

        for rank_value, tied_methods in sorted(rank_groups.items()):
            if len(tied_methods) < 2:
                continue

            for method_a, method_b in itertools.combinations(sorted(tied_methods), 2):
                record = {
                    "replicate": replicate,
                    "rank": rank_value,
                    "method_a": method_a,
                    "method_b": method_b,
                }

                row_a = row_lookup.get((replicate, method_a))
                row_b = row_lookup.get((replicate, method_b))

                if row_a is None or row_b is None:
                    record["action"] = "missing_cmpnets_dataframe_row"
                    comparison_records.append(record)
                    continue

                record["false_positive_a"] = row_a["false_positive"]
                record["false_negative_a"] = row_a["false_negative"]
                record["network_file_a"] = row_a["network_file"]

                record["false_positive_b"] = row_b["false_positive"]
                record["false_negative_b"] = row_b["false_negative"]
                record["network_file_b"] = row_b["network_file"]

                if find(parent, method_a) == find(parent, method_b):
                    record["action"] = "already_same_network"
                    record["same_network"] = True
                    comparison_records.append(record)
                    continue

                fp_equal = is_close(row_a["false_positive"], row_b["false_positive"])
                fn_equal = is_close(row_a["false_negative"], row_b["false_negative"])

                if not (fp_equal and fn_equal):
                    record["action"] = "fp_fn_not_equal"
                    record["same_network"] = False
                    comparison_records.append(record)
                    continue

                rep_a = find(parent, method_a)
                rep_b = find(parent, method_b)

                true_method, computed_method = sorted([rep_a, rep_b])
                record["true_method"] = true_method
                record["computed_method"] = computed_method

                pair_key = (true_method, computed_method)

                if pair_key in checked_representative_pairs:
                    record["action"] = "representative_already_checked"
                    record["same_network"] = checked_representative_pairs[pair_key]
                    comparison_records.append(record)
                    continue

                true_row = row_lookup.get((replicate, true_method))
                computed_row = row_lookup.get((replicate, computed_method))

                if true_row is None or computed_row is None:
                    record["action"] = "missing_representative_row"
                    comparison_records.append(record)
                    continue

                true_fp_fn_equal = (
                    is_close(true_row["false_positive"], computed_row["false_positive"])
                    and is_close(true_row["false_negative"], computed_row["false_negative"])
                )

                if not true_fp_fn_equal:
                    record["action"] = "representative_fp_fn_not_equal"
                    record["same_network"] = False
                    checked_representative_pairs[pair_key] = False
                    comparison_records.append(record)
                    continue

                pair_dir = os.path.join(
                    output_path,
                    replicate,
                    f"{slugify(true_method)}__vs__{slugify(computed_method)}",
                )
                os.makedirs(pair_dir, exist_ok=True)

                nexus_file, cmpnets_output_file = create_nexus_file(
                    pair_dir,
                    true_row["network_file"],
                    computed_row["network_file"],
                    mode,
                )

                log_file = os.path.join(pair_dir, "phylonet.log")

                record["pair_dir"] = pair_dir
                record["nexus_file"] = nexus_file
                record["cmpnets_output_file"] = cmpnets_output_file
                record["log_file"] = log_file

                if dry_run:
                    record["action"] = "dry_run"
                    record["same_network"] = None
                    comparison_records.append(record)
                    continue

                command = f'java -jar "{phylonet_path}" "{nexus_file}" > "{log_file}" 2>&1'
                os.system(command)

                record["action"] = "ran_phylonet"

                if not os.path.exists(cmpnets_output_file):
                    record["action"] = "missing_cmpnets_output"
                    record["same_network"] = None
                    comparison_records.append(record)
                    continue

                fn, fp, fa = get_fn_fp_fa(cmpnets_output_file)

                record["cmpnets_false_negative"] = fn
                record["cmpnets_false_positive"] = fp
                record["cmpnets_false_average"] = fa

                same_network = is_close(fn, 0.0) and is_close(fp, 0.0)

                record["same_network"] = same_network
                checked_representative_pairs[pair_key] = same_network

                if same_network:
                    union(parent, true_method, computed_method)

                comparison_records.append(record)

        group_row = {"replicate": replicate}
        root_labels = {}

        for method in methods:
            root = find(parent, method)
            if root not in root_labels:
                root_labels[root] = len(root_labels) + 1
            group_row[method] = root_labels[root]

        group_records.append(group_row)

    comparisons = pd.DataFrame(comparison_records)
    groups = pd.DataFrame(group_records)

    if not groups.empty:
        groups = groups[["replicate"] + methods]

    comparisons_path = os.path.join(output_path, "same_rank_network_comparisons.csv")
    groups_path = os.path.join(output_path, "same_network_groups.csv")

    comparisons.to_csv(comparisons_path, index=False)
    groups.to_csv(groups_path, index=False)

    print(f"Wrote {comparisons_path}")
    print(f"Wrote {groups_path}")

    if not groups.empty:
        print("\nSame-network groups (same number in a row = same network):")
        print(groups.to_string(index=False))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script to compare same-rank inferred networks using Cmpnets",
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

    input_csv = config["input_csv"]
    ranks_csv = config["ranks_csv"]
    output_path = config["output_path"]
    phylonet_path = config["phylonet_path"]
    mode = config["mode"]
    dry_run = config["dry_run"]

    compare_same_rank_networks(input_csv,ranks_csv,output_path,phylonet_path,mode,dry_run)
