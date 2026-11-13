#!/usr/bin/env python3
"""Re-run gene-tree inference only for replicates that are missing it.

Scans ``<input>/subset_*/*/*/`` for replicate directories without an
``iqtree_500.nwk`` and runs ``infer_gene_trees.py`` on each. Because the
inference is resumable (genes with an existing treefile are skipped), this
finishes the partial replicates without redoing completed ones::

    python3 redo_missing.py --input data --threads 32

With ``--list`` it only prints the missing replicate paths, one per line,
for use as a SLURM job array::

    python3 redo_missing.py --input data --list
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

GENE_TREE_NAME = "iqtree_500.nwk"


def find_missing(input_dir: Path, name: str) -> list[Path]:
    """Replicate directories under ``input_dir`` without a ``name`` file."""
    return [d for d in sorted(input_dir.glob("subset_*/*/*/"))
            if not (d / name).exists()]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Infer gene trees only for replicates missing them.")
    parser.add_argument("--input", type=Path, default=Path("data"),
                        help="directory containing the subset_* folders")
    parser.add_argument("--threads", type=int, default=0,
                        help="number of concurrent IQ-TREE processes "
                             "(default: all cores)")
    parser.add_argument("--name", default=GENE_TREE_NAME,
                        help=f"gene-tree file name (default: {GENE_TREE_NAME})")
    parser.add_argument("--list", action="store_true",
                        help="only print the missing replicate paths, one "
                             "per line (for a SLURM job array)")
    args = parser.parse_args(argv)

    input_dir = args.input.resolve()
    todo = find_missing(input_dir, args.name)
    if not todo:
        print("nothing to do: every replicate has its gene trees")
        return

    if args.list:
        for d in todo:
            print(d)
        return

    print(f"{len(todo)} replicate(s) missing {args.name}:")
    for d in todo:
        print(f"  {d.relative_to(input_dir)}")

    cmd = [sys.executable,
           str(Path(__file__).resolve().parent / "infer_gene_trees.py"),
           "--name", args.name]
    if args.threads > 0:
        cmd += ["--threads", str(args.threads)]

    failed = 0
    for i, rep in enumerate(todo, 1):
        rel = rep.relative_to(input_dir)
        print(f"[{i}/{len(todo)}] {rel}", flush=True)
        result = subprocess.run(cmd + ["--input", str(rep)])
        if result.returncode != 0:
            failed += 1
            print(f"[{i}/{len(todo)}] {rel}: FAILED", file=sys.stderr)

    if failed:
        sys.exit(f"{failed} of {len(todo)} replicates failed")


if __name__ == "__main__":
    main()
