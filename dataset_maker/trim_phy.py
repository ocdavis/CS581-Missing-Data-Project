#!/usr/bin/env python3
"""Trim PHYLIP (.phy) sequence datasets to fixed sequence lengths.

Reads every ``.phy`` file under an input directory, parses it into a simple
in-memory structure (header + list of ``(name, sequence)`` rows), truncates
each sequence to a target length, and writes the result to an output
directory while preserving the relative file structure.

A config file (TOML) supplies a list of lengths; one full copy of the dataset
is produced per length::

    output_dir / f"len_{length}" / <relative path of the input file>

Usage::

    python3 trim_phy.py --config example_config.toml
    python3 trim_phy.py --input data/camus-dataset/n50 --output data/trimmed \
        --lengths 100 250 500 --workers 8

Files are processed in parallel with a process pool (``--workers``, default:
all cores). Each file is parsed once and written once per length. A progress
bar is shown when tqdm is installed; otherwise a plain counter is printed.

The parsing / trimming / writing steps are separate functions so they can be
reused or overridden from other scripts.
"""

from __future__ import annotations

import argparse
import os
import sys
import tomllib
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

try:
    from tqdm import tqdm
except ImportError:  # graceful fallback: simple counter instead of a bar
    tqdm = None

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
                 lengths: Sequence[int]) -> list[Path]:
    """Parse one file and write one trimmed copy per length.

    The file is parsed only once regardless of how many lengths are requested.
    Returns the list of destination paths.
    """
    rel = src.relative_to(input_dir)
    data = parse_phylip(src)
    dsts = []
    for length in lengths:
        dst = output_dir / f"len_{length}" / rel
        write_phylip(trim_sequences(data, length), dst)
        dsts.append(dst)
    return dsts


def _progress(iterable, total: int, desc: str):
    """Wrap ``iterable`` with tqdm if available, else yield items unchanged."""
    if tqdm is not None:
        return tqdm(iterable, total=total, desc=desc, unit="file")
    return iterable


def _print_fallback_progress(done: int, total: int, desc: str) -> None:
    if tqdm is None and (done % 25 == 0 or done == total):
        print(f"\r{desc}: {done}/{total} files", end="" if done < total else "\n",
              file=sys.stderr)


def run(input_dir: Path, output_dir: Path, lengths: Sequence[int],
        workers: int | None = None) -> None:
    files = list(find_phy_files(input_dir))
    if not files:
        print(f"warning: no .phy files found under {input_dir}", file=sys.stderr)
        return

    if workers is None:
        workers = os.cpu_count() or 1
    workers = max(1, min(workers, len(files)))

    desc = f"trimming ({len(lengths)} lengths)"
    if workers == 1:
        for done, src in enumerate(_progress(files, len(files), desc), 1):
            process_file(src, input_dir, output_dir, lengths)
            _print_fallback_progress(done, len(files), desc)
        return

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(process_file, src, input_dir, output_dir, lengths)
                   for src in files]
        for done, future in enumerate(
                _progress(as_completed(futures), len(files), desc), 1):
            future.result()  # propagate worker exceptions
            _print_fallback_progress(done, len(files), desc)


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
    parser.add_argument("--workers", type=int,
                        help="parallel worker processes (default: all cores)")
    args = parser.parse_args(argv)

    if args.config:
        cfg = load_config(args.config)
        input_dir = Path(cfg["input_dir"])
        output_dir = Path(cfg["output_dir"])
        lengths = [int(x) for x in cfg["lengths"]]
        workers = cfg.get("workers")
    elif args.input and args.output and args.lengths:
        input_dir, output_dir, lengths = args.input, args.output, args.lengths
        workers = None
    else:
        parser.error("provide either --config or all of --input/--output/--lengths")

    if args.workers is not None:
        workers = args.workers

    run(input_dir, output_dir, lengths, workers)


if __name__ == "__main__":
    main()
