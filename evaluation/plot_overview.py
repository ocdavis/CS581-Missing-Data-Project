#!/usr/bin/env python3
"""One-page overview of every experiment: FN and FP box plots per method.

One row per experiment (sequence length, taxa missing from every gene, and
taxa missing from a share of the genes), FN rate on the left and FP rate on
the right, all on the same scale. Reads the CSV files written by
compare_tob.py:

    python3 evaluation/plot_overview.py results/base_tree_errors.csv \\
        results/tob_qmc_errors.csv results/ectoblob_errors.csv \\
        --output results/overview.png

Values above the top of the scale (trees with no internal branches have a
false negative rate of 100%) are drawn as triangles at the top edge. Needs
matplotlib.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_errors import GRID, INK, METHODS, MUTED

# (row title, x-axis label, [(condition folder, tick label), ...])
ROWS = [
    ("Gene tree estimation error", "Sequence length (sites)", [
        ("miss_0_len_500", "500"), ("miss_0_len_200", "200"),
        ("miss_0_len_100", "100"), ("miss_0_len_50_files_100", "50"),
        ("miss_0_len_25_files_100", "25")]),
    ("Taxa missing from every gene", "Taxa missing", [
        ("miss_0_len_500", "0%"), ("miss_20_len_500", "20%"),
        ("miss_50_len_500", "50%"), ("miss_80_len_500_files_100", "80%")]),
    ("50% of taxa missing from some genes", "Genes with 50% of taxa missing", [
        ("miss_0_len_500", "0%"), ("miss_50_len_500_files_50", "50%"),
        ("miss_50_len_500_files_80", "80%"), ("miss_50_len_500_files_95", "95%"),
        ("miss_50_len_500_files_99", "99%"), ("miss_50_len_500", "100%")]),
    ("80% of taxa missing from some genes", "Genes with 80% of taxa missing", [
        ("miss_0_len_500", "0%"), ("miss_80_len_500_files_50", "50%"),
        ("miss_80_len_500_files_80", "80%"), ("miss_80_len_500_files_95", "95%"),
        ("miss_80_len_500_files_99", "99%"),
        ("miss_80_len_500_files_100", "100%")]),
]
PANELS = [("fn_rate", "False negative rate (%)"),
          ("fp_rate", "False positive rate (%)")]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("csv", nargs="+", type=Path,
                        help="CSV file(s) written by compare_tob.py")
    parser.add_argument("--output", type=Path, default=Path("overview.png"))
    parser.add_argument("--top", type=float, default=70,
                        help="top of the rate scale in percent (default: 70)")
    args = parser.parse_args()

    # values[column][method label][condition] -> rates in percent
    values = {col: defaultdict(lambda: defaultdict(list)) for col, _ in PANELS}
    for path in args.csv:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                label = next((lab for pat, lab, _ in METHODS
                              if fnmatch.fnmatch(row["method"], pat)), None)
                if label:
                    for col, _ in PANELS:
                        values[col][label][row["condition"]].append(
                            float(row[col]) * 100)
    methods = [(lab, colour) for _, lab, colour in METHODS
               if lab in values["fn_rate"]]

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 11,
                         "axes.edgecolor": MUTED, "text.color": INK,
                         "axes.labelcolor": MUTED, "xtick.color": MUTED,
                         "ytick.color": MUTED})
    fig, axes = plt.subplots(len(ROWS), 2, figsize=(15, 4.1 * len(ROWS)),
                             sharey=True, squeeze=False)
    step = 0.8 / len(methods)
    off_scale = False
    for r, (title, xlabel, conditions) in enumerate(ROWS):
        for c, (col, ylabel) in enumerate(PANELS):
            ax = axes[r][c]
            for k, (label, colour) in enumerate(methods):
                xs, data, high = [], [], []
                for i, (cond, _) in enumerate(conditions):
                    vals = values[col][label].get(cond)
                    if vals:
                        x = i - 0.4 + (k + 0.5) * step
                        xs.append(x)
                        data.append([v for v in vals if v <= args.top])
                        high += [x] * sum(v > args.top for v in vals)
                if not xs:
                    continue
                box = ax.boxplot(
                    data, positions=xs, widths=step * 0.8, patch_artist=True,
                    manage_ticks=False,
                    boxprops={"facecolor": colour, "edgecolor": INK,
                              "linewidth": 0.8},
                    medianprops={"color": INK, "linewidth": 1.4},
                    whiskerprops={"color": INK, "linewidth": 0.8},
                    capprops={"color": INK, "linewidth": 0.8},
                    flierprops={"marker": "o", "markersize": 3,
                                "markerfacecolor": colour,
                                "markeredgecolor": INK,
                                "markeredgewidth": 0.4})
                box["boxes"][0].set_label(label)
                if high:
                    off_scale = True
                    ax.plot(high, [args.top * 0.985] * len(high), "^",
                            color=colour, markeredgecolor=INK,
                            markeredgewidth=0.4, markersize=7, clip_on=False)
            ax.set_xticks(range(len(conditions)), [t for _, t in conditions])
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.tick_params(labelleft=True, length=0)
            ax.yaxis.grid(True, color=GRID, linewidth=0.8)
            ax.set_axisbelow(True)
            ax.spines[["top", "right", "left"]].set_visible(False)
        axes[r][0].set_title(title, loc="left", fontsize=13, color=INK,
                             fontweight="bold", pad=10)
    axes[0][0].set_ylim(0, args.top)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, ncols=len(methods), loc="upper right",
               frameon=False, handlelength=1.2, fontsize=12)
    fig.suptitle("FN and FP rates of the estimated trees of blobs, CAMUS n50",
                 x=0.01, ha="left", fontsize=15, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97), h_pad=2.0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=150, facecolor="white")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
