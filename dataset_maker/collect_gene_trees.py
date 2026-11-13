#!/usr/bin/env python3
"""Collect the inferred gene trees from every subset into one directory.

Every ``iqtree_500.nwk`` under ``data/subset_*/<dataset>/<rep>/`` is copied
to ``<output_dir>/<dataset>/<rep>/iqtree_500.nwk``, so the result has the
same layout as ``data/trimmed`` but with all replicates together::

    python3 collect_gene_trees.py --input data --output data/gene_trees

Files are copied with ``shutil.copy2`` (metadata preserved) and are skipped
when the destination already exists, so the script is resumable. The exit
status is non-zero if any subset replicate is missing its gene-tree file.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

GENE_TREE_NAME = "iqtree_500.nwk"


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Collect iqtree_500.nwk from every subset into one "
                    "directory with the trimmed layout.")
    parser.add_argument("--input", type=Path, default=Path("data"),
                        help="directory containing the subset_* folders")
    parser.add_argument("--output", type=Path, default=Path("data/gene_trees"),
                        help="output directory (default: data/gene_trees)")
    parser.add_argument("--name", default=GENE_TREE_NAME,
                        help=f"gene-tree file name (default: {GENE_TREE_NAME})")
    args = parser.parse_args(argv)

    input_dir = args.input.resolve()
    output_dir = args.output.resolve()

    # <subset>/<dataset>/<rep>/<name>, replicates sorted numerically.
    sources = sorted(input_dir.glob(f"subset_*/*/*/{args.name}"),
                     key=lambda p: (p.parent.parent.name, p.parent.name))
    if not sources:
        sys.exit(f"no {args.name} files under {input_dir}/subset_*")

    copied = skipped = 0
    missing: list[Path] = []
    for src in sources:
        dataset, rep = src.parent.parent.name, src.parent.name
        dest = output_dir / dataset / rep / args.name
        if dest.exists():
            skipped += 1
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        copied += 1

    # Report replicates that exist in the subsets but have no gene trees.
    for rep_dir in sorted(input_dir.glob("subset_*/*/*/")):
        if not (rep_dir / args.name).exists():
            missing.append(rep_dir.relative_to(input_dir))

    print(f"done: {copied} copied, {skipped} already present, "
          f"{len(missing)} missing")
    if missing:
        for rel in missing:
            print(f"  missing: {rel}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
