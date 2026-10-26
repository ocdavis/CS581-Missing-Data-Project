import re
import tomllib
import os
import pandas as pd

def remove_branch_lengths(net):
    x = re.sub(r":[0-9.:]*(?=[,)])", "", net)
    return re.sub(r"(?<=\))[0-9.:]+(?=[,;)])", "", x)


def get_mpl_score(path):
    with open(path, "r") as file:
        for line_number, line in enumerate(file, start=1):
            if line_number == 5:
                return float(line.split(":", 1)[1].strip())

def get_fn_fp(path):
    with open(path, "r") as file:
        for line_number, line in enumerate(file, start=1):
            if line_number == 2:
                info = (line.split(":", 1)[1].strip()).split()
                return float(info[0]), float(info[1])

def get_fn_fp_fa(path):
    with open(path, "r") as file:
        for line_number, line in enumerate(file, start=1):
            if line_number == 2:
                info = (line.split(":", 1)[1].strip()).split()
                return float(info[0]), float(info[1]), float(info[2])

def load_config(path, required_keys):
    with open(path, "rb") as config_file:
        config = tomllib.load(config_file)

    missing_keys = required_keys - config.keys()

    if missing_keys:
        missing = ", ".join(sorted(missing_keys))
        raise ValueError(
            f"Configuration file is missing required keys: {missing}"
        )

    return config

def merge_records_with_stats(df_records, df_stats, how='inner'):
    """Merge accuracy metrics with performance metrics."""
    merged = df_records.merge(
        df_stats,
        left_on=['file_path', 'replicate'],
        right_on=['output_file', 'replicate'],
        how=how
    )
    # Drop the redundant output_file column
    merged = merged.drop('file_path', axis=1)
    return merged

def cmpnets_consolidate_records(input_path, file_names, method_labels):
    records = []

    for root, dirs, files in os.walk(input_path):
        replicate_name = os.path.relpath(root, input_path)
        for file_name in files:
            if file_name in file_names:
                file_path = os.path.join(root, file_name)
                fn, fp, fa = get_fn_fp_fa(file_path)
                records.append({
                    "method": file_name,
                    "file_path": file_path,
                    "false_negative": fn,
                    "false_positive": fp,
                    "false_average": fa,
                    "replicate": replicate_name,
                })

    df = pd.DataFrame(records)
    df["method_label"] = df["method"].map(method_labels)

    return df
