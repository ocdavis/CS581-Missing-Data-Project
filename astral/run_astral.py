#!/usr/bin/env python3
"""Infer an ASTRAL species tree for every gene-tree file under a directory.

For each <dir>/.../iqtree_500.nwk it runs ASTRAL with its default options

    java -jar astral.5.7.8.jar -i <genes> -o <out>

and writes the species tree next to the gene trees:

    data/gene_trees/miss_0_len_500/00/iqtree_500.nwk
      -> data/gene_trees/miss_0_len_500/00/astral.nwk

With --output the trees are written under another directory instead, keeping
the path relative to the input directory:

    --output data/astral_trees
      -> data/astral_trees/miss_0_len_500/00/astral.nwk

ASTRAL's messages go to <out stem>.log (astral.log). Existing trees are
skipped, so an interrupted run can simply be restarted.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_JAR = PROJECT / "repos" / "ASTRAL" / "Astral" / "astral.5.7.8.jar"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", nargs="?", type=Path,
                        default=Path("data/gene_trees"),
                        help="directory searched recursively for gene-tree files "
                             "(default: data/gene_trees)")
    parser.add_argument("--output", type=Path,
                        help="write the species trees under this directory, "
                             "keeping each file's path relative to the input "
                             "directory (default: next to the gene trees)")
    parser.add_argument("--pattern", default="iqtree_500.nwk",
                        help="gene-tree file name to look for "
                             "(default: iqtree_500.nwk)")
    parser.add_argument("--out-name", default="astral.nwk",
                        help="species-tree file name (default: astral.nwk)")
    parser.add_argument("--jar", type=Path, default=DEFAULT_JAR,
                        help=f"ASTRAL jar (default: {DEFAULT_JAR})")
    parser.add_argument("--overwrite", action="store_true",
                        help="rerun ASTRAL even where a species tree exists")
    return parser.parse_args()


def out_path(genes: Path, args: argparse.Namespace) -> Path:
    if args.output is None:
        return genes.parent / args.out_name
    return args.output / genes.relative_to(args.input).parent / args.out_name


def run_one(genes: Path, out: Path, jar: Path) -> bool:
    """Run ASTRAL for one gene-tree file; True if a tree was written."""
    log = out.with_suffix(".log")
    tmp = out.with_name(out.name + ".tmp")
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["java", "-jar", str(jar), "-i", str(genes), "-o", str(tmp)]
    with open(log, "w") as logfile:
        ok = subprocess.run(cmd, stdout=logfile, stderr=logfile).returncode == 0
    if ok and tmp.exists() and tmp.stat().st_size > 0:
        tmp.replace(out)
        return True
    tmp.unlink(missing_ok=True)
    return False


def main() -> int:
    args = parse_args()
    if not args.input.is_dir():
        print(f"error: {args.input} is not a directory", file=sys.stderr)
        return 2
    if not args.jar.is_file():
        print(f"error: {args.jar} not found", file=sys.stderr)
        return 2

    genes = sorted(args.input.rglob(args.pattern))
    if not genes:
        print(f"error: no {args.pattern} files under {args.input}",
              file=sys.stderr)
        return 2

    print(f"=== ASTRAL for {len(genes)} gene-tree file(s) under {args.input}, "
          f"jar {args.jar} ===")
    n_written = n_skipped = 0
    failed = []
    for i, gene_file in enumerate(genes, 1):
        out = out_path(gene_file, args)
        if out.exists() and out.stat().st_size > 0 and not args.overwrite:
            n_skipped += 1
            print(f"[{i}/{len(genes)}] skip   {out}")
            continue
        print(f"[{i}/{len(genes)}] run    {gene_file}")
        if run_one(gene_file, out, args.jar):
            n_written += 1
            print(f"[{i}/{len(genes)}] wrote  {out}")
        else:
            failed.append(gene_file)
            print(f"[{i}/{len(genes)}] FAILED {gene_file} "
                  f"(see {out.with_suffix('.log')})", file=sys.stderr)

    print(f"=== done: {n_written} written, {n_skipped} already present, "
          f"{len(failed)} failed ===")
    for gene_file in failed:
        print(f"FAILED {gene_file}", file=sys.stderr)
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
