#!/usr/bin/env python3
"""Polytomies of the estimated trees of blobs next to the true ones.

One row per experiment (the same rows as plot_overview.py). Left: mean
number of polytomies per tree, with 95% confidence intervals. Right: box
plots of the largest polytomy degree (3 = no polytomy). The true trees of
blobs are drawn as their own series. Reads the CSV files written by
polytomy_stats.py:

    python3 evaluation/plot_polytomies.py results/tob_qmc_polytomies.csv \\
        results/ectoblob_polytomies.csv --output results/polytomies.jpg

Needs matplotlib.
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
from plot_overview import ROWS

TRUE = ("True tree of blobs", "#52514e")


def mean_ci(values: list[float]) -> tuple[float, float]:
    n = len(values)
    mean = sum(values) / n
    if n < 2:
        return mean, 0.0
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    return mean, 1.96 * (var / n) ** 0.5


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("csv", nargs="+", type=Path,
                        help="CSV file(s) written by polytomy_stats.py")
    parser.add_argument("--output", type=Path, default=Path("polytomies.jpg"))
    args = parser.parse_args()

    # values[statistic][series label][condition] -> one value per replicate
    values = {s: defaultdict(lambda: defaultdict(dict))
              for s in ("polytomies", "max_degree")}
    for path in args.csv:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                label = next((lab for pat, lab, _ in METHODS
                              if fnmatch.fnmatch(row["method"], pat)), None)
                if not label:
                    continue
                cond, rep = row["condition"], row["replicate"]
                for s in values:
                    values[s][label][cond][rep] = float(row[f"est_{s}"])
                    values[s][TRUE[0]][cond][rep] = float(row[f"ref_{s}"])
    series = [TRUE] + [(lab, colour) for _, lab, colour in METHODS
                       if lab in values["polytomies"]]

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 11,
                         "axes.edgecolor": MUTED, "text.color": INK,
                         "axes.labelcolor": MUTED, "xtick.color": MUTED,
                         "ytick.color": MUTED})
    fig, axes = plt.subplots(len(ROWS), 2, figsize=(15, 4.1 * len(ROWS)),
                             sharey="col", squeeze=False)
    step = 0.8 / len(series)
    for r, (title, xlabel, conditions) in enumerate(ROWS):
        for c, stat in enumerate(("polytomies", "max_degree")):
            ax = axes[r][c]
            for k, (label, colour) in enumerate(series):
                xs, data = [], []
                for i, (cond, _) in enumerate(conditions):
                    vals = list(values[stat][label].get(cond, {}).values())
                    if vals:
                        xs.append(i - 0.4 + (k + 0.5) * step)
                        data.append(vals)
                if not xs:
                    continue
                if stat == "polytomies":
                    means, cis = zip(*(mean_ci(v) for v in data))
                    ax.bar(xs, means, width=step * 0.9, color=colour,
                           label=label, yerr=cis, capsize=2,
                           error_kw={"ecolor": INK, "elinewidth": 0.8})
                else:
                    box = ax.boxplot(
                        data, positions=xs, widths=step * 0.8,
                        patch_artist=True, manage_ticks=False,
                        boxprops={"facecolor": colour, "edgecolor": INK,
                                  "linewidth": 0.8},
                        medianprops={"color": "white" if label == TRUE[0]
                                     else INK, "linewidth": 1.4},
                        whiskerprops={"color": INK, "linewidth": 0.8},
                        capprops={"color": INK, "linewidth": 0.8},
                        flierprops={"marker": "o", "markersize": 3,
                                    "markerfacecolor": colour,
                                    "markeredgecolor": INK,
                                    "markeredgewidth": 0.4})
                    box["boxes"][0].set_label(label)
            ax.set_xticks(range(len(conditions)), [t for _, t in conditions])
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Polytomies per tree (mean)" if stat == "polytomies"
                          else "Largest polytomy degree")
            ax.tick_params(labelleft=True, length=0)
            ax.yaxis.grid(True, color=GRID, linewidth=0.8)
            ax.set_axisbelow(True)
            ax.spines[["top", "right", "left"]].set_visible(False)
        axes[r][0].set_title(title, loc="left", fontsize=13, color=INK,
                             fontweight="bold", pad=10)
    axes[0][0].set_ylim(bottom=0)
    axes[0][1].set_ylim(bottom=0)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, ncols=len(series), loc="upper right",
               frameon=False, handlelength=1.2, fontsize=12)
    fig.suptitle("Polytomies of the estimated and true trees of blobs, "
                 "CAMUS n50", x=0.01, ha="left", fontsize=15, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97), h_pad=2.0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=150, facecolor="white")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
