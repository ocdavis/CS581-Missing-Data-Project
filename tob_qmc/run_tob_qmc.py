#!/usr/bin/env python3
"""Run TOB-QMC (``tree-qmc --blob``) on every gene-tree file in a directory.

For each file matching ``gene_tree_glob`` under ``input_dir`` this runs::

    tree-qmc --blob --iter_limit_blob 0 --alpha 1e-7 --beta 0.95 --root OUT \
        -i <gene trees> -o <output_dir>/<relative path>/<name>.tob.tre

and saves tree-qmc's messages next to the tree as ``<name>.tob.log``.
Files whose output tree already exists are skipped.

Usage::

    python3 run_tob_qmc.py --config example_config.toml
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tomllib
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", type=Path, required=True,
                        help="TOML config file (see example_config.toml)")
    args = parser.parse_args()

    with args.config.open("rb") as fh:
        cfg = tomllib.load(fh)
    tree_qmc = cfg.get("tree_qmc", "tree-qmc")
    input_dir = Path(cfg["input_dir"])
    output_dir = Path(cfg["output_dir"])
    alpha = cfg.get("alpha", 1e-7)
    beta = cfg.get("beta", 0.95)
    outgroup = cfg.get("outgroup", "OUT")
    exhaustive = cfg.get("exhaustive", True)

    files = sorted(input_dir.rglob(cfg["gene_tree_glob"]))
    if not files:
        sys.exit(f"no files matching {cfg['gene_tree_glob']!r} under {input_dir}")

    failed = 0
    for i, genes in enumerate(files, 1):
        rel = genes.relative_to(input_dir)
        out = output_dir / rel.parent / f"{rel.stem}.tob.tre"
        log = out.with_suffix(".log")
        if out.exists():
            print(f"[{i}/{len(files)}] {rel}: already done, skipping")
            continue

        cmd = [tree_qmc, "--blob", "--alpha", f"{alpha:g}", "--beta", f"{beta:g}",
               "--root", outgroup, "-i", str(genes), "-o", str(out), "--override"]
        if exhaustive:
            cmd += ["--iter_limit_blob", "0"]

        out.parent.mkdir(parents=True, exist_ok=True)
        with log.open("w") as fh:
            result = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            failed += 1
            out.unlink(missing_ok=True)
            print(f"[{i}/{len(files)}] {rel}: FAILED, see {log}")
        else:
            print(f"[{i}/{len(files)}] {rel}: wrote {out}")

    if failed:
        sys.exit(f"{failed} of {len(files)} files failed")


if __name__ == "__main__":
    main()
