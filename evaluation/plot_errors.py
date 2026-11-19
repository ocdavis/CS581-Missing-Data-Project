#!/usr/bin/env python3
"""Box plots of FN and FP per condition from compare_tob.py CSV files.

One group of boxes per condition, one box per method (the ``method`` column,
i.e. the file name of the estimated trees), in two panels: false negatives
and false positives. Each box covers the replicates of one condition: the
line is the median, the box the quartiles, the whiskers reach the furthest
value within 1.5 box lengths, and dots are values beyond that.

    python3 evaluation/plot_errors.py results/base_tree_errors.csv \\
        results/tob_qmc_errors.csv --output results/fn_fp.png

Any number of CSV files can be given, so ECToBlob's is simply added to the
list. Rates (percent) are plotted by default; ``--counts`` plots numbers of
branches instead. Needs matplotlib.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# File name pattern -> (label, colour). Listed in drawing order; a method not
# listed here keeps its file name and takes the next spare colour.
METHODS = [
    ("tree-qmc*.nwk", "TREE-QMC base tree", "#a3a29b"),
    ("tob-qmc*.nwk", "TOB-QMC", "#2a78d6"),
    ("*ectoblob*", "ECToBlob", "#eb6834"),
]
SPARE_COLOURS = ["#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


CONDITION = re.compile(r"miss_(\d+)_len_(\d+)(?:_files_(\d+))?")


def condition_key(name: str) -> tuple:
    """Sort miss_<m>_len_<l>[_files_<f>] conditions: original first, then
    shorter sequences, then more missing taxa (in fewer genes first). Other
    names sort alphabetically."""
    m = CONDITION.fullmatch(name)
    if not m:
        return (1, name)
    miss, length, files = int(m.group(1)), int(m.group(2)), int(m.group(3) or 100)
    return (0, miss, -length, files)


def condition_label(name: str, longest: int, one_length: bool) -> str:
    m = CONDITION.fullmatch(name)
    if not m:
        return name or "all"
    miss, length = int(m.group(1)), int(m.group(2))
    if miss == 0:
        # when only the missing taxa vary, name the complete data by that
        return "0% missing" if one_length else f"{length} bp"
    label = f"{miss}% missing"
    if m.group(3) and m.group(3) != "100":
        label += f"\nin {m.group(3)}% of genes"
    if length != longest:
        label = f"{length} bp\n" + label
    return label


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("csv", nargs="+", type=Path,
                        help="CSV file(s) written by compare_tob.py")
    parser.add_argument("--output", type=Path, default=Path("fn_fp.png"),
                        help="image file; the extension picks the format "
                             "(.png, .pdf, .svg) (default: fn_fp.png)")
    parser.add_argument("--counts", action="store_true",
                        help="plot numbers of branches instead of rates")
    parser.add_argument("--conditions", nargs="+", metavar="NAME",
                        help="only plot these conditions (as named in the "
                             "CSV, e.g. miss_0_len_500 miss_0_len_100)")
    parser.add_argument("--xlabel",
                        help="label under the x axis of each panel")
    parser.add_argument("--title", default="",
                        help="title above the figure")
    args = parser.parse_args()

    # values[panel][method][condition] -> list over replicates
    panels = ("fn", "fp") if args.counts else ("fn_rate", "fp_rate")
    scale = 1 if args.counts else 100
    values = {p: defaultdict(lambda: defaultdict(list)) for p in panels}
    for path in args.csv:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if args.conditions and row["condition"] not in args.conditions:
                    continue
                for p in panels:
                    values[p][row["method"]][row["condition"]].append(
                        float(row[p]) * scale)

    found = list(values[panels[0]])
    if not found:
        raise SystemExit("no rows to plot")
    methods, spare = [], iter(SPARE_COLOURS)
    for pattern, label, colour in METHODS:
        for name in [m for m in found if fnmatch.fnmatch(m, pattern)]:
            if name not in [m[0] for m in methods]:
                # several files of one method (e.g. other alphas) keep their names
                first = not any(m[1] == label for m in methods)
                methods.append((name, label if first else name,
                                colour if first else next(spare)))
    for name in found:
        if name not in [m[0] for m in methods]:
            methods.append((name, name, next(spare)))

    conditions = sorted({c for m in found for c in values[panels[0]][m]},
                        key=condition_key)
    lengths = [int(m.group(2)) for c in conditions
               if (m := CONDITION.fullmatch(c))]
    longest = max(lengths) if lengths else 0
    one_length = len(conditions) > 1 and len(set(lengths)) == 1
    labels = [condition_label(c, longest, one_length) for c in conditions]
    xlabel = args.xlabel
    parsed = [CONDITION.fullmatch(c) for c in conditions]
    if one_length and all(parsed):
        # (percent of taxa missing, percent of genes affected) per condition;
        # no _files_ suffix means every gene
        pairs = [(int(m.group(1)), int(m.group(3) or 100) if int(m.group(1)) else 0)
                 for m in parsed]
        missing = {miss for miss, _ in pairs if miss}
        shares = {share for miss, share in pairs if miss}
        if shares == {100}:
            # only the percentage of missing taxa varies: put it on the axis
            labels = [f"{miss}%" for miss, _ in pairs]
            if xlabel is None:
                xlabel = "Taxa missing from every gene"
        elif len(missing) == 1:
            # only the share of incomplete genes varies: put that on the axis
            labels = [f"{share}%" for _, share in pairs]
            if xlabel is None:
                xlabel = f"Genes with {missing.pop()}% of taxa missing"

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 12,
                         "axes.edgecolor": MUTED, "text.color": INK,
                         "axes.labelcolor": MUTED, "xtick.color": MUTED,
                         "ytick.color": MUTED})
    # few conditions: the two panels side by side; many: one above the other
    if len(conditions) <= 3:
        fig, axes = plt.subplots(1, 2, figsize=(max(9.0, 3.0 * len(conditions)), 4.2),
                                 sharey=True)
    else:
        fig, axes = plt.subplots(2, 1, figsize=(2.3 * len(conditions), 7.4),
                                 sharex=True, sharey=True)
    if args.counts:
        titles = {"fn": "False negative branches per tree",
                  "fp": "False positive branches per tree"}
    else:
        titles = {"fn": "False negative rate (%)",
                  "fp": "False positive rate (%)"}
    step = 0.8 / len(methods)
    for ax, panel in zip(axes, panels):
        for k, (name, label, colour) in enumerate(methods):
            xs, data = [], []
            for i, cond in enumerate(conditions):
                vals = values[panel][name].get(cond)
                if vals:
                    xs.append(i - 0.4 + (k + 0.5) * step)
                    data.append(vals)
            box = ax.boxplot(
                data, positions=xs, widths=step * 0.8, patch_artist=True,
                manage_ticks=False,
                boxprops={"facecolor": colour, "edgecolor": INK, "linewidth": 0.8},
                medianprops={"color": INK, "linewidth": 1.4},
                whiskerprops={"color": INK, "linewidth": 0.8},
                capprops={"color": INK, "linewidth": 0.8},
                flierprops={"marker": "o", "markersize": 3,
                            "markerfacecolor": colour,
                            "markeredgecolor": INK, "markeredgewidth": 0.4})
            box["boxes"][0].set_label(label)
        ax.tick_params(labelleft=True)
        ax.set_ylabel(titles[panel[:2]])
        ax.set_xticks(range(len(conditions)), labels)
        if xlabel and (len(conditions) <= 3 or ax is axes[-1]):
            ax.set_xlabel(xlabel)
        ax.yaxis.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(length=0)
    # one scale for both panels, so FN and FP can be compared directly
    top = max(v for p in panels for m in values[p].values()
              for c in conditions for v in m.get(c, []))
    axes[0].set_ylim(0, top * 1.05 if top else 1)
    fig.legend(*axes[0].get_legend_handles_labels(), ncols=len(methods),
               loc="upper right", frameon=False, handlelength=1.2)
    if args.title:
        fig.suptitle(args.title, x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200, facecolor="white")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
