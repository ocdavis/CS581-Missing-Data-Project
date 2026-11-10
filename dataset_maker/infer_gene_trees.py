#!/usr/bin/env python3
"""Infer gene trees with IQ-TREE for every replicate of a dataset directory.

For each ``seqs-500`` folder found under the input directory (e.g.
``data/trimmed/len_100/00/seqs-500``), every ``.phy`` alignment inside is
given to IQ-TREE and the resulting ML trees are collected into a single
multi-tree Newick file, one tree per line, written next to the ``seqs-500``
folder (or under ``--output``, mirroring the directory layout)::

    python3 infer_gene_trees.py --input data/trimmed/len_100
    python3 infer_gene_trees.py --input data/camus-dataset/n50 --threads 16

The output file is named ``iqtree_500.nwk`` by default, matching the gene-tree
files distributed with the CAMUS datasets, regardless of the actual
alignment length.

IQ-TREE is run one alignment at a time (``-T 1``) with ``--threads``
processes running concurrently; for the many small alignments in these
datasets this is faster than giving IQ-TREE multiple threads per run.
With ``--model MFP``, the best-fit model is selected on the
first gene of each replicate and that model is then reused for every other
gene of the replicate; any other ``--model`` value is used as-is for all
genes. A ``--seed`` makes runs reproducible.

Already-computed genes are skipped: a gene is only (re-)run when its
``<stem>.treefile`` does not exist yet, so an interrupted run can simply be
restarted. Pass ``--overwrite`` to redo everything. Intermediate IQ-TREE
files (``.log``, ``.iqtree``, ``.ckp.gz``, ...) are kept in an ``iqtree-work``
folder next to each ``seqs-500`` folder (or under ``--workdir DIR``,
mirroring the input layout); ``--clean`` deletes the work folder of each
replicate that finished successfully.

A replicate only gets its output file when every alignment in it produced a
tree; otherwise the replicate is reported as failed and left without one.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tomllib
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

DEFAULT_NAME = "iqtree_500.nwk"
DEFAULT_MODEL = "GTR+G"
DEFAULT_SEED = 581
DEFAULT_WORKDIR_NAME = "iqtree-work"


def find_replicates(input_dir: Path) -> list[Path]:
    """Return every ``seqs-500`` directory under ``input_dir``, sorted."""
    return sorted(p for p in input_dir.rglob("seqs-500") if p.is_dir())


def treefile_for(prefix: Path) -> Path:
    """Path of the IQ-TREE ML tree for the run with prefix ``prefix``."""
    return prefix.with_suffix(".treefile")


def best_fit_model(prefix: Path) -> str | None:
    """Model chosen by MFP in the IQ-TREE report ``prefix.iqtree``, if any."""
    report = prefix.with_suffix(".iqtree")
    if not report.exists():
        return None
    for line in report.read_text().splitlines():
        if "Best-fit model according to BIC" in line:
            return line.split(":")[-1].strip()
    return None


def run_iqtree(phy: Path, prefix: Path, model: str, seed: int,
               iqtree_bin: str) -> Path:
    """Run IQ-TREE on one alignment and return the path of its ML tree.

    All IQ-TREE output (``.treefile``, ``.iqtree``, ``.log``, ...) is written
    under ``prefix``. Raises ``RuntimeError`` with the tail of IQ-TREE's
    output on failure.
    """
    treefile = treefile_for(prefix)
    # -redo: ignore any leftover checkpoint from an interrupted earlier run
    # (we only get here when the treefile itself is missing).
    cmd = [iqtree_bin, "-s", str(phy), "-pre", str(prefix),
           "-m", model, "-T", "1", "-seed", str(seed), "-redo"]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True)
    if proc.returncode != 0 or not treefile.exists():
        tail = "\n".join(proc.stdout.splitlines()[-15:])
        raise RuntimeError(f"IQ-TREE failed on {phy} (exit {proc.returncode}):\n{tail}")
    return treefile


def process_replicate(seqs_dir: Path, out_file: Path, model: str, seed: int,
                      threads: int, overwrite: bool, iqtree_bin: str,
                      workdir: Path, clean: bool) -> bool:
    """Infer trees for every alignment in ``seqs_dir`` and write ``out_file``.

    IQ-TREE output prefixes live under ``workdir`` (by default an
    ``iqtree-work`` folder next to ``seqs_dir``). With ``clean`` the workdir
    is deleted after a successful run. Returns True on success, False if any
    alignment failed or has no tree.
    """
    phys = sorted(seqs_dir.glob("*.phy"))
    if not phys:
        print(f"WARNING: no .phy files in {seqs_dir}", file=sys.stderr)
        return False

    def prefix_for(phy: Path) -> Path:
        return workdir / phy.stem

    for prefix in (prefix_for(p) for p in phys):
        prefix.parent.mkdir(parents=True, exist_ok=True)

    todo = [p for p in phys if overwrite or not treefile_for(prefix_for(p)).exists()]
    failures: list[str] = []

    # With MFP, select the best-fit model on the first gene only and reuse
    # it for every other gene of the replicate; a fixed model is used as-is.
    effective = model
    if model.upper() == "MFP":
        first = phys[0]
        effective = best_fit_model(prefix_for(first))
        if effective is None:
            try:
                run_iqtree(first, prefix_for(first), "MFP", seed, iqtree_bin)
            except RuntimeError as exc:
                failures.append(str(exc))
            effective = best_fit_model(prefix_for(first))
            if effective is None:
                print(f"ERROR: could not determine best-fit model from "
                      f"{first.name} in {seqs_dir}", file=sys.stderr)
                return False
        print(f"{out_file.parent.name}: best-fit model from {first.name}: "
              f"{effective}", file=sys.stderr)
        todo = [p for p in todo if p is not first]

    if todo:
        is_tty = sys.stderr.isatty()
        pbar = (tqdm(total=len(todo), desc=str(out_file.parent.name), unit="gene")
                if tqdm is not None and is_tty else None)
        done = 0
        with ThreadPoolExecutor(max_workers=threads) as pool:
            futures = {pool.submit(run_iqtree, p, prefix_for(p), effective,
                                   seed, iqtree_bin): p
                       for p in todo}
            pending = set(futures)
            while pending:
                finished, pending = wait(pending, return_when=FIRST_COMPLETED)
                for fut in finished:
                    phy = futures[fut]
                    try:
                        fut.result()
                    except RuntimeError as exc:
                        failures.append(str(exc))
                    done += 1
                    if pbar is not None:
                        pbar.update(1)
                    elif is_tty:
                        print(f"\r{out_file.parent.name}: {done}/{len(todo)}",
                              end="", file=sys.stderr, flush=True)
        if pbar is not None:
            pbar.close()
        elif is_tty:
            print(file=sys.stderr)

    trees = []
    missing = []
    for phy in phys:
        tf = treefile_for(prefix_for(phy))
        if tf.exists():
            trees.append(tf.read_text().strip())
        else:
            missing.append(phy.name)
    if missing or failures:
        for name in missing:
            print(f"ERROR: no tree for {name} in {seqs_dir}", file=sys.stderr)
        for msg in failures:
            print(msg, file=sys.stderr)
        print(f"ERROR: {out_file} not written "
              f"({len(missing) + len(failures)} of {len(phys)} genes failed)",
              file=sys.stderr)
        return False

    out_file.write_text("\n".join(trees) + "\n")
    if clean:
        shutil.rmtree(workdir, ignore_errors=True)
    return True


def load_config(config_path: Path) -> dict:
    with config_path.open("rb") as fh:
        cfg = tomllib.load(fh)
    if "input_dir" not in cfg:
        raise ValueError(f"config {config_path} is missing key 'input_dir'")
    return cfg


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Infer gene trees with IQ-TREE for every seqs-500 folder "
                    "under an input directory.")
    parser.add_argument("--config", type=Path,
                        help="TOML config file (input_dir, output_dir, model, "
                             "seed, threads, name)")
    parser.add_argument("--input", type=Path,
                        help="input directory (overrides config)")
    parser.add_argument("--output", type=Path,
                        help="output directory; by default the gene-tree file "
                             "is written next to each seqs-500 folder")
    parser.add_argument("--model", default=DEFAULT_MODEL,
                        help=f"IQ-TREE model (default: {DEFAULT_MODEL})")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help=f"IQ-TREE seed (default: {DEFAULT_SEED})")
    parser.add_argument("--threads", type=int, default=0,
                        help="number of concurrent IQ-TREE processes "
                             "(default: all cores)")
    parser.add_argument("--name", default=DEFAULT_NAME,
                        help=f"output file name (default: {DEFAULT_NAME})")
    parser.add_argument("--iqtree", default="iqtree",
                        help="iqtree binary to call (default: iqtree on PATH)")
    parser.add_argument("--workdir", type=Path,
                        help="directory for IQ-TREE's intermediate files "
                             "(.treefile, .log, ...), mirroring the input "
                             "layout; by default they go into an "
                             f"'{DEFAULT_WORKDIR_NAME}' folder next to each "
                             "seqs-500 folder")
    parser.add_argument("--clean", action="store_true",
                        help="delete each replicate's work directory after "
                             "it finishes successfully")
    parser.add_argument("--overwrite", action="store_true",
                        help="re-run genes whose treefile already exists")
    args = parser.parse_args(argv)

    if args.config:
        cfg = load_config(args.config)
        if args.input is None:
            args.input = Path(cfg["input_dir"])
        if args.output is None and "output_dir" in cfg:
            args.output = Path(cfg["output_dir"])
        if args.workdir is None and "workdir" in cfg:
            args.workdir = Path(cfg["workdir"])
        if "clean" in cfg:
            args.clean = cfg["clean"]
        if args.model == DEFAULT_MODEL and "model" in cfg:
            args.model = cfg["model"]
        if args.seed == DEFAULT_SEED and "seed" in cfg:
            args.seed = cfg["seed"]
        if args.threads == 0 and "threads" in cfg:
            args.threads = cfg["threads"]
        if args.name == DEFAULT_NAME and "name" in cfg:
            args.name = cfg["name"]
    if args.input is None:
        parser.error("provide either --config or --input")

    threads = args.threads if args.threads > 0 else (os.cpu_count() or 1)
    input_dir = args.input.resolve()
    global_workdir = args.workdir.resolve() if args.workdir else None
    replicates = find_replicates(input_dir)
    if not replicates:
        parser.error(f"no seqs-500 folders found under {input_dir}")

    print(f"{len(replicates)} replicate(s), model {args.model}, "
          f"seed {args.seed}, {threads} concurrent job(s)"
          + (", cleaning work directories" if args.clean else ""))

    ok = 0
    failed: list[Path] = []
    for seqs_dir in replicates:
        rel = seqs_dir.parent.relative_to(input_dir)
        out_file = ((args.output.resolve() / rel) if args.output
                    else seqs_dir.parent) / args.name
        out_file.parent.mkdir(parents=True, exist_ok=True)
        workdir = ((global_workdir / rel) if global_workdir
                   else seqs_dir.parent / DEFAULT_WORKDIR_NAME)
        if process_replicate(seqs_dir, out_file, args.model, args.seed,
                             threads, args.overwrite, args.iqtree, workdir,
                             args.clean):
            ok += 1
        else:
            failed.append(rel)

    print(f"done: {ok} written, {len(failed)} failed")
    if failed:
        for rel in failed:
            print(f"  failed: {rel}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
