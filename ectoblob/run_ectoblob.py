#!/usr/bin/env python3
"""Run ECToBlob (MSCquartets R package) on every gene-tree file in a directory.

For each file matching ``gene_tree_glob`` under ``input_dir`` this looks up
the starting tree (a resolution of the true tree of blobs, e.g. a TREE-QMC
species tree) matching ``start_tree_glob`` in the same directory and runs::

    Rscript ectoblob.R --genes <gene trees> --tree <start tree> \
        --out <output_dir>/<relative path>/<name>.ectoblob.tre [options]

saving R's messages next to the tree as ``<name>.ectoblob.log``. Files whose
output tree already exists are skipped.

Usage::

    python3 run_ectoblob.py --config example_config.toml
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
    rscript = cfg.get("rscript", "Rscript")
    r_script = Path(cfg.get("r_script", Path(__file__).parent / "ectoblob.R"))
    input_dir = Path(cfg["input_dir"])
    gene_tree_glob = cfg["gene_tree_glob"]
    start_tree_glob = cfg["start_tree_glob"]
    output_dir = Path(cfg["output_dir"])
    alpha = cfg.get("alpha", 0.05)
    beta = cfg.get("beta", 0.8)
    qtype = cfg.get("qtype", "bi")
    test_correction = cfg.get("test_correction", "Bon")
    index = cfg.get("index", "early")
    outgroup = cfg.get("outgroup", "OUT")
    save_table = cfg.get("save_table", False)

    files = sorted(input_dir.rglob(gene_tree_glob))
    if not files:
        sys.exit(f"no files matching {gene_tree_glob!r} under {input_dir}")

    failed = 0
    for i, genes in enumerate(files, 1):
        rel = genes.relative_to(input_dir)
        out = output_dir / rel.parent / f"{rel.stem}.ectoblob.tre"
        log = out.with_suffix(".log")
        if out.exists():
            print(f"[{i}/{len(files)}] {rel}: already done, skipping")
            continue

        # The starting tree must be a resolution of the true tree of blobs
        # for this replicate (e.g. a TREE-QMC species tree); exactly one
        # match of start_tree_glob is expected next to the gene trees.
        starts = sorted(genes.parent.glob(start_tree_glob))
        if len(starts) != 1:
            failed += 1
            print(f"[{i}/{len(files)}] {rel}: FAILED, {len(starts)} files "
                  f"match {start_tree_glob!r} in {genes.parent} (need 1)")
            continue

        cmd = [rscript, str(r_script),
               "--genes", str(genes), "--tree", str(starts[0]),
               "--out", str(out),
               "--alpha", f"{alpha:g}", "--beta", f"{beta:g}",
               "--qtype", qtype, "--correction", test_correction,
               "--index", index]
        if outgroup:
            cmd += ["--outgroup", outgroup]
        if save_table:
            cmd += ["--save-table"]

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
