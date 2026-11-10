#!/usr/bin/env python3
"""Trim PHYLIP (.phy) sequence datasets to fixed sequence lengths.

Reads every ``.phy`` file under an input directory, parses it into a simple
in-memory structure (header + list of ``(name, sequence)`` rows), truncates
each sequence to a target length, and writes the result to an output
directory while preserving the relative file structure.

A config file (TOML) supplies a list of lengths; one full copy of the dataset
is produced per length::

    output_dir / f"len_{length}" / <relative path of the input file>

It can also supply a list of missing-taxa percentages. The two lists are then
crossed: for every (percentage, length) pair, one copy of the dataset is
produced in which every file has that percentage of its ingroup taxa
(everything except the outgroup) removed and its sequences truncated to that
length::

    output_dir / f"miss_{percent}" / f"len_{length}" / <relative path>

An independent random subset of taxa is drawn for each file, and the same
subset is used at every length. Include 0 in the percentages to get copies
with all taxa, and the original sequence length in the lengths to get
untrimmed copies; ``miss_0/len_<original length>`` is then a copy of the
input alignments, so every condition lives in the same directory layout.

Usage::

    python3 trim_phy.py --config example_config.toml
    python3 trim_phy.py --input data/camus-dataset/n50 --output data/trimmed \
        --lengths 100 250 500 --workers 8
    python3 trim_phy.py --input data/camus-dataset/n50 --output data/trimmed \
        --lengths 100 250 500 --missing 20 50 --outgroup OUT --seed 581

Files are processed in parallel with a process pool (``--workers``, default:
all cores). Each file is parsed once and written once per length. Files are
streamed from the directory walk with only a bounded window of tasks in
flight, so memory stays flat regardless of dataset size. Progress is shown
as a live bar on a terminal, or as periodic percentage lines when the output
is piped or redirected.

The parsing / trimming / writing steps are separate functions so they can be
reused or overridden from other scripts.
"""

from __future__ import annotations

import argparse
import os
import random
import sys
import tomllib
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from typing import Iterable, Iterator, Sequence

try:
    from tqdm import tqdm
except ImportError:  # graceful fallback: simple counter instead of a bar
    tqdm = None

DEFAULT_OUTGROUP = "OUT"  # outgroup taxon name in the CAMUS dataset
DEFAULT_SEED = 0

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class PhylipFile:
    """In-memory representation of a PHYLIP interleaved/sequential file."""

    n_taxa: int
    n_sites: int
    rows: list[tuple[str, str]]  # (taxon name, sequence)

    def __post_init__(self) -> None:
        if len(self.rows) != self.n_taxa:
            raise ValueError(
                f"header says {self.n_taxa} taxa but parsed {len(self.rows)} rows"
            )


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def parse_phylip(path: Path) -> PhylipFile:
    """Parse a PHYLIP file (sequential or interleaved) into a PhylipFile.

    Handles both layouts:
      * sequential: each taxon's full sequence follows its name (possibly
        wrapped over continuation lines that start with whitespace);
      * interleaved: the taxon name is repeated at the start of every block.
    """
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError(f"{path}: empty file")

    header = lines[0].split()
    if len(header) != 2 or not all(part.isdigit() for part in header):
        raise ValueError(f"{path}: malformed header line {lines[0]!r}")
    n_taxa, n_sites = int(header[0]), int(header[1])

    rows: dict[str, str] = {}
    order: list[str] = []
    last_name: str | None = None

    for line in lines[1:]:
        if not line.strip():
            continue
        if line[0].isspace():
            # Continuation line (no taxon name) -> append to previous taxon.
            if last_name is None:
                raise ValueError(f"{path}: continuation line before any taxon")
            rows[last_name] += line.strip()
        else:
            name, _, seq = line.partition(" ")
            name = name.strip()
            seq = seq.strip()
            if name in rows:
                # Interleaved block: name repeats, append the new chunk.
                rows[name] += seq
            else:
                rows[name] = seq
                order.append(name)
                last_name = name

    return PhylipFile(n_taxa=n_taxa, n_sites=n_sites,
                      rows=[(name, rows[name]) for name in order])


# ---------------------------------------------------------------------------
# Transformation (swap this out or compose it to add your own improvements)
# ---------------------------------------------------------------------------


def trim_sequences(data: PhylipFile, length: int) -> PhylipFile:
    """Return a new PhylipFile with every sequence truncated to ``length``.

    Sequences shorter than ``length`` are left as-is (no padding).
    """
    new_rows = [(name, seq[:length]) for name, seq in data.rows]
    new_sites = min(length, data.n_sites)
    return PhylipFile(n_taxa=data.n_taxa, n_sites=new_sites, rows=new_rows)


def n_taxa_to_remove(n_ingroup: int, percent: float) -> int:
    """Number of ingroup taxa to drop for a given missing percentage."""
    return round(n_ingroup * percent / 100)


def remove_random_taxa(data: PhylipFile, percent: float, rng: random.Random,
                       outgroup: str | None = DEFAULT_OUTGROUP) -> PhylipFile:
    """Return a new PhylipFile with ``percent`` % of the ingroup taxa removed.

    The taxa to drop are sampled uniformly without replacement using ``rng``.
    The outgroup is never removed and is not counted when computing how many
    taxa to drop (e.g. 20% of 50 ingroup taxa + outgroup removes 10 taxa and
    leaves 41 rows). Remaining rows keep their original order. Pass
    ``outgroup=None`` to treat every taxon as removable.
    """
    if not 0 <= percent < 100:
        raise ValueError(f"missing percentage must be in [0, 100): {percent}")
    names = [name for name, _ in data.rows]
    if outgroup is not None and outgroup not in names:
        raise ValueError(f"outgroup {outgroup!r} not found among taxa")
    ingroup = [name for name in names if name != outgroup]
    dropped = set(rng.sample(ingroup, n_taxa_to_remove(len(ingroup), percent)))
    new_rows = [row for row in data.rows if row[0] not in dropped]
    return PhylipFile(n_taxa=len(new_rows), n_sites=data.n_sites,
                      rows=new_rows)


def file_rng(seed: int, rel: Path, percent: float) -> random.Random:
    """Deterministic RNG for one (file, missing percentage) pair.

    Seeding from the file's relative path (rather than sharing one RNG) makes
    the result independent of worker count and processing order, so reruns
    with the same seed reproduce the same dataset.
    """
    return random.Random(f"{seed}|{rel.as_posix()}|{percent:g}")


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------


def write_phylip(data: PhylipFile, path: Path) -> None:
    """Write a PhylipFile in sequential (one line per taxon) PHYLIP format."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(f"{data.n_taxa}  {data.n_sites}\n")
        for name, seq in data.rows:
            fh.write(f"{name:<10} {seq}\n")


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def find_phy_files(input_dir: Path) -> Iterable[Path]:
    """Yield every .phy file under ``input_dir`` (recursively)."""
    yield from sorted(input_dir.rglob("*.phy"))


def process_file(src: Path, input_dir: Path, output_dir: Path,
                 lengths: Sequence[int],
                 missing: Sequence[float] = (),
                 outgroup: str | None = DEFAULT_OUTGROUP,
                 seed: int = DEFAULT_SEED) -> list[Path]:
    """Parse one file and write one copy per (missing percentage, length).

    The file is parsed only once regardless of how many variants are
    requested. Without ``missing``, copies go to ``len_<L>/``. With
    ``missing``, every percentage is crossed with every length and copies go
    to ``miss_<P>/len_<L>/``; the taxa removed for a given percentage are the
    same at every length. Returns the list of destination paths.
    """
    rel = src.relative_to(input_dir)
    data = parse_phylip(src)
    dsts = []
    if not missing:
        for length in lengths:
            dst = output_dir / f"len_{length}" / rel
            write_phylip(trim_sequences(data, length), dst)
            dsts.append(dst)
        return dsts
    for percent in missing:
        try:
            reduced = remove_random_taxa(
                data, percent, file_rng(seed, rel, percent), outgroup)
        except ValueError as exc:
            raise ValueError(f"{src}: {exc}") from exc
        for length in lengths:
            dst = output_dir / f"miss_{percent:g}" / f"len_{length}" / rel
            write_phylip(trim_sequences(reduced, length), dst)
            dsts.append(dst)
    return dsts


def run(input_dir: Path, output_dir: Path, lengths: Sequence[int],
        workers: int | None = None,
        missing: Sequence[float] = (),
        outgroup: str | None = DEFAULT_OUTGROUP,
        seed: int = DEFAULT_SEED) -> None:
    """Stream files from ``input_dir`` and process them in parallel.

    Files are never all held in memory: the directory walk is consumed
    lazily and only ``workers * 4`` tasks are in flight at any time. A quick
    metadata-only pre-walk counts the files up front so the progress bar can
    show a percentage.

    Progress reporting adapts to the output: a live tqdm bar when stderr is
    a terminal, and periodic one-line-per-update progress messages when the
    output is piped or redirected (where ``\\r``-based bars cannot render).
    """
    total = sum(1 for _ in find_phy_files(input_dir))
    if total == 0:
        print(f"warning: no .phy files found under {input_dir}", file=sys.stderr)
        return

    if workers is None:
        workers = os.cpu_count() or 1
    workers = max(1, min(workers, total))

    desc = (f"processing {total} files ({len(lengths)} lengths, "
            f"{len(missing)} missing levels)")
    is_tty = sys.stderr.isatty()
    pbar = (tqdm(total=total, desc=desc, unit="file")
            if tqdm is not None and is_tty else None)
    done_count = 0
    last_reported = 0
    # When piped, report every ~2% so logs show incremental progress
    # without flooding them.
    report_step = max(1, total // 50)

    def advance(n: int) -> None:
        nonlocal done_count, last_reported
        done_count += n
        if pbar is not None:
            pbar.update(n)
        elif is_tty:
            # No tqdm installed: simple in-place counter on a terminal.
            print(f"\r{desc}: {done_count}/{total} files",
                  end="" if done_count < total else "\n",
                  file=sys.stderr, flush=True)
        elif done_count - last_reported >= report_step or done_count == total:
            # Piped/redirected: newline-delimited updates that survive logs.
            pct = 100 * done_count / total
            print(f"{desc}: {done_count}/{total} files ({pct:.0f}%)",
                  file=sys.stderr, flush=True)
            last_reported = done_count

    try:
        if workers == 1:
            for src in find_phy_files(input_dir):
                process_file(src, input_dir, output_dir, lengths,
                             missing, outgroup, seed)
                advance(1)
            return

        # Bounded window of in-flight futures; the walk is consumed lazily.
        max_in_flight = workers * 4
        file_iter: Iterator[Path] = iter(find_phy_files(input_dir))
        with ProcessPoolExecutor(max_workers=workers) as pool:
            pending = {
                pool.submit(process_file, src, input_dir, output_dir, lengths,
                            missing, outgroup, seed)
                for src in islice(file_iter, max_in_flight)
            }
            while pending:
                finished, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in finished:
                    future.result()  # propagate worker exceptions
                advance(len(finished))
                for src in islice(file_iter, len(finished)):
                    pending.add(
                        pool.submit(process_file, src, input_dir, output_dir,
                                    lengths, missing, outgroup, seed))
    finally:
        if pbar is not None:
            pbar.close()


def load_config(config_path: Path) -> dict:
    with config_path.open("rb") as fh:
        cfg = tomllib.load(fh)
    for key in ("input_dir", "output_dir", "lengths"):
        if key not in cfg:
            raise ValueError(f"config {config_path} is missing key {key!r}")
    return cfg


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", type=Path,
                        help="TOML config file (see example_config.toml)")
    parser.add_argument("--input", type=Path,
                        help="input directory (overrides config)")
    parser.add_argument("--output", type=Path,
                        help="output directory (overrides config)")
    parser.add_argument("--lengths", type=int, nargs="+",
                        help="sequence lengths (overrides config)")
    parser.add_argument("--missing", type=float, nargs="+",
                        help="percentages of ingroup taxa to remove per file "
                             "(overrides config)")
    parser.add_argument("--outgroup",
                        help="taxon that is never removed (overrides config)")
    parser.add_argument("--seed", type=int,
                        help="random seed for taxon removal (overrides config)")
    parser.add_argument("--workers", type=int,
                        help="parallel worker processes (default: all cores)")
    args = parser.parse_args(argv)

    if args.config:
        cfg = load_config(args.config)
        input_dir = Path(cfg["input_dir"])
        output_dir = Path(cfg["output_dir"])
        lengths = [int(x) for x in cfg["lengths"]]
        missing = [float(x) for x in cfg.get("missing_percents", [])]
        outgroup = cfg.get("outgroup", DEFAULT_OUTGROUP)
        seed = cfg.get("seed", DEFAULT_SEED)
        workers = cfg.get("workers")
    elif args.input and args.output and args.lengths:
        input_dir, output_dir, lengths = args.input, args.output, args.lengths
        missing, outgroup, seed = [], DEFAULT_OUTGROUP, DEFAULT_SEED
        workers = None
    else:
        parser.error("provide either --config or all of --input/--output/--lengths")

    if args.missing is not None:
        missing = args.missing
    if args.outgroup is not None:
        outgroup = args.outgroup
    if args.seed is not None:
        seed = args.seed
    if args.workers is not None:
        workers = args.workers

    run(input_dir, output_dir, lengths, workers, missing, outgroup, seed)


if __name__ == "__main__":
    main()
