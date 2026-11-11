#!/usr/bin/env python3
"""Generate the nine trimmed/missing-taxa datasets and infer gene trees
for a subset of replicates.

Three steps:

1. run ``trim_phy.py`` with ``example_config.toml`` to produce the nine
   datasets ``data/trimmed/miss_<P>_len_<L>/`` (P in {0, 20, 50}, L in
   {100, 200, 500}), unless ``data/trimmed`` already exists — so of two
   simultaneously started runs only the first generates the datasets;
2. copy the given replicates of the five datasets of interest into
   ``data/subset_<range>/``, preserving the layout (skipped if the subset
   directory already exists):
     miss_0_len_100, miss_0_len_200, miss_0_len_500  (0% missing, all lengths)
     miss_20_len_500, miss_50_len_500                (missing, full length)
3. run ``infer_gene_trees.py`` on the subset (``iqtree_500.nwk`` written
   next to each ``seqs-500`` folder; finished genes are skipped on rerun).

Usage::

    python3 run_subset.py 00-04 --threads 32
    python3 run_subset.py 05-09
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
DATASET_MAKER = PROJECT / "dataset_maker"

# The five datasets of interest: 0% missing at every length, plus the
# missing-taxa datasets at the original (unshortened) length.
DATASETS = [
    "miss_0_len_100",
    "miss_0_len_200",
    "miss_0_len_500",
    "miss_20_len_500",
    "miss_50_len_500",
]


def parse_range(text: str) -> list[str]:
    """Expand a replicate range like ``00-04`` into ``00 01 02 03 04``."""
    start, sep, end = text.partition("-")
    if not sep:
        return [start]
    return [f"{i:02d}" for i in range(int(start), int(end) + 1)]


def run(cmd: list[str]) -> None:
    print(f"=== running: {' '.join(cmd)} ===", flush=True)
    subprocess.run(cmd, check=True)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Generate the nine datasets and infer gene trees for a "
                    "subset of replicates.")
    parser.add_argument("range",
                        help="replicate range, e.g. 00-04 or 05-09 "
                             "(a single replicate like 07 also works)")
    parser.add_argument("--threads", type=int, default=0,
                        help="number of concurrent IQ-TREE processes "
                             "(default: all cores)")
    parser.add_argument("--trimmed-dir", type=Path, default=PROJECT / "data" / "trimmed",
                        help="directory the nine datasets are generated in")
    parser.add_argument("--subset-root", type=Path, default=PROJECT / "data",
                        help="directory the subset is created under")
    args = parser.parse_args(argv)

    replicates = parse_range(args.range)
    subset = args.subset_root / f"subset_{args.range.replace('-', '_')}"

    # --- step 1: the nine datasets (only if not generated yet) --------------
    if not args.trimmed_dir.exists():
        # Create the directory first so a simultaneously started run skips
        # trim_phy instead of racing us on the same output files.
        args.trimmed_dir.mkdir(parents=True)
        run([sys.executable, str(DATASET_MAKER / "trim_phy.py"),
             "--config", str(DATASET_MAKER / "example_config.toml")])
    else:
        print(f"=== {args.trimmed_dir} exists, skipping trim_phy ===")

    # --- step 2: assemble the subset (only if not assembled yet) ------------
    if not subset.exists():
        print(f"=== copying replicates {' '.join(replicates)} into {subset} ===")
        for ds in DATASETS:
            for rep in replicates:
                dest = subset / ds / rep
                dest.mkdir(parents=True)
                shutil.copytree(args.trimmed_dir / ds / rep / "seqs-500",
                                dest / "seqs-500")
    else:
        print(f"=== {subset} exists, skipping copy ===")

    # --- step 3: gene trees --------------------------------------------------
    cmd = [sys.executable, str(DATASET_MAKER / "infer_gene_trees.py"),
           "--input", str(subset), "--clean"]
    if args.threads > 0:
        cmd += ["--threads", str(args.threads)]
    run(cmd)


if __name__ == "__main__":
    main()
