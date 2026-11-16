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

    output_dir / f"miss_{percent}_len_{length}" / <relative path>

An independent random subset of taxa is drawn for each file, and the same
subset is used at every length. Include 0 in the percentages to get copies
with all taxa, and the original sequence length in the lengths to get
untrimmed copies; ``miss_0_len_<original length>`` is then a copy of the
input alignments, so every condition lives in the same directory layout.

A third list, ``file_percents``, controls how many of the files are affected.
It is crossed with both lists, so for every (file percentage, missing
percentage, length) triple one copy of the dataset is produced in which only
that percentage of the files have taxa removed; the remaining files keep all
their taxa and are merely truncated to the requested length::

    output_dir / f"miss_{percent}_files_{files}_len_{length}" / <relative path>

For example ``--lengths 500 --missing 50 --file-percents 50`` writes
``miss_50_files_50_len_500``: half the alignments lose half their ingroup
taxa, the other half stay complete, and every sequence keeps its original
length. A file percentage of 100 reduces every file, which is the default and
keeps the plain ``miss_<P>_len_<L>`` names used by the existing datasets.

The number of files picked is exact, not sampled: the files are grouped by
their parent directory (one group per ``seqs-500`` replicate) and, within each
group, the ``round(n * F / 100)`` files with the lowest seeded draw are
chosen. Every replicate is therefore affected exactly ``F`` % of the time
whenever that lands on a whole gene, and the counts are identical whatever
``--workers`` is. The number of taxa removed is likewise exact: ``round`` of
``n_ingroup * P / 100``, drawn uniformly without replacement from the ingroup.

Usage::

    python3 trim_phy.py --config example_config.toml
    python3 trim_phy.py --input data/camus-dataset/n50 --output data/trimmed \
        --lengths 100 250 500 --workers 8
    python3 trim_phy.py --input data/camus-dataset/n50 --output data/trimmed \
        --lengths 100 250 500 --missing 20 50 --outgroup OUT --seed 581
    python3 trim_phy.py --input data/camus-dataset/n50 --output data/trimmed \
        --lengths 500 --missing 50 --file-percents 50

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
from typing import Iterable, Iterator, Mapping, Sequence

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


def _round_count(n: int, percent: float) -> int:
    """``n`` scaled by ``percent``, rounded to the nearest whole item.

    Used both for the number of taxa to drop and the number of files to
    select, so the two behave the same way.
    """
    return round(n * percent / 100)


def n_taxa_to_remove(n_ingroup: int, percent: float) -> int:
    """Number of ingroup taxa to drop for a given missing percentage."""
    return _round_count(n_ingroup, percent)


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


def file_draw(seed: int, rel: Path) -> float:
    """Uniform draw in [0, 1) that identifies one file.

    Seeded from the file's relative path with its own salt, so it is
    independent of the taxon draw (:func:`file_rng`) and of the number of
    workers or the processing order. The draw deliberately does not depend on
    any percentage, which is what lets :func:`selection_cutoffs` rank files
    and keeps the conditions nested.
    """
    return random.Random(f"{seed}|{rel.as_posix()}|filemiss").random()


def selection_cutoffs(files: Iterable[Path], input_dir: Path, seed: int,
                      file_percent: float) -> dict[Path, float]:
    """Per-directory thresholds that select an *exact* number of files.

    Files are grouped by their parent directory — one group per ``seqs-500``
    folder in the CAMUS layout, or a single group for a flat input directory.
    Within each group the files are ranked by :func:`file_draw` and the
    ``round(n * file_percent / 100)`` lowest are selected, so every group
    gets exactly the requested count whenever the percentage lands on a whole
    file and the nearest whole number otherwise. Ranking each directory
    separately also means every replicate of a dataset is affected equally
    often, instead of the count fluctuating from replicate to replicate.

    Returns ``{parent directory: threshold}``; a file is selected when its
    draw is ``<=`` the threshold of its own directory. Because the ranking
    uses the same draw at every percentage, the set selected at ``F`` is a
    subset of the set selected at any ``F' > F``.
    """
    if not 0 < file_percent <= 100:
        raise ValueError(f"file percentage must be in (0, 100]: {file_percent}")
    groups: dict[Path, list[tuple[float, Path]]] = {}
    for src in files:
        rel = src.relative_to(input_dir)
        groups.setdefault(rel.parent, []).append((file_draw(seed, rel), rel))
    cutoffs: dict[Path, float] = {}
    for parent, items in groups.items():
        items.sort()
        n_select = _round_count(len(items), file_percent)
        if n_select <= 0:
            cutoffs[parent] = float("-inf")     # select nothing
        elif n_select >= len(items):
            cutoffs[parent] = float("inf")      # select everything
        else:
            cutoffs[parent] = items[n_select - 1][0]
    return cutoffs


def file_selected(seed: int, rel: Path, file_percent: float,
                  cutoff: float | None = None) -> bool:
    """Deterministic decision: does ``rel`` get any taxa removed?

    ``file_percent`` is the percentage of the files that should be reduced;
    100 selects every file. With ``cutoff`` — a threshold from
    :func:`selection_cutoffs` — the file is selected when its draw is at or
    below it, which makes the number of selected files exact. Without a
    cutoff the draw is compared against ``file_percent`` directly, which is
    correct in expectation but only approximately ``file_percent`` % of the
    files; that path exists for callers that process files one at a time and
    cannot see the whole file set.

    Either way the draw does not depend on ``file_percent``, so the
    selections for different percentages are *nested*: every file reduced at
    ``F`` is also reduced at any ``F' > F``, and the affected files are the
    same at every missing percentage and length. Comparing two conditions
    therefore isolates how many genes are affected rather than also changing
    which genes they are. Like :func:`file_rng`, the choice is independent of
    worker count and processing order.
    """
    if not 0 < file_percent <= 100:
        raise ValueError(f"file percentage must be in (0, 100]: {file_percent}")
    if file_percent == 100:
        return True
    draw = file_draw(seed, rel)
    if cutoff is not None:
        return draw <= cutoff
    return draw < file_percent / 100


def condition_dir(missing_percent: float, length: int,
                  file_percent: float = 100) -> str:
    """Directory name for one (file %, missing %, length) condition.

    The ``files_<F>`` component is left out when every file is reduced
    (``F == 100``) so the datasets written before this option existed keep
    their ``miss_<P>_len_<L>`` names.
    """
    if file_percent >= 100:
        return f"miss_{missing_percent:g}_len_{length}"
    return f"miss_{missing_percent:g}_files_{file_percent:g}_len_{length}"


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
                 seed: int = DEFAULT_SEED,
                 file_percents: Sequence[float] = (),
                 cutoffs: Mapping[float, Mapping[Path, float]] | None = None
                 ) -> list[Path]:
    """Parse one file and write one copy per condition.

    The file is parsed only once regardless of how many variants are
    requested. Without ``missing``, copies go to ``len_<L>/``. With
    ``missing``, every file percentage is crossed with every missing
    percentage and every length, and copies go to
    ``miss_<P>_files_<F>_len_<L>/`` (plain ``miss_<P>_len_<L>/`` when ``F``
    is 100, the default). A file percentage of ``F`` removes taxa only from
    the files picked by :func:`file_selected`; the other files keep all their
    taxa and are merely truncated to the requested length. The taxa removed
    for a given percentage are the same at every length. Returns the list of
    destination paths.

    ``cutoffs`` maps a file percentage to the threshold table from
    :func:`selection_cutoffs`. Passing it makes the number of selected files
    exact; leaving it out falls back to the per-file comparison, which is
    only approximately ``F`` % of the files.
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
    for file_percent in (tuple(file_percents) or (100.0,)):
        thresholds = None if cutoffs is None else cutoffs.get(file_percent)
        selected = file_selected(
            seed, rel, file_percent,
            None if thresholds is None else thresholds[rel.parent])
        for percent in missing:
            if selected:
                try:
                    reduced = remove_random_taxa(
                        data, percent, file_rng(seed, rel, percent), outgroup)
                except ValueError as exc:
                    raise ValueError(f"{src}: {exc}") from exc
            else:
                # Not chosen for reduction: keep every taxon, trim only.
                reduced = data
            for length in lengths:
                dst = (output_dir
                       / condition_dir(percent, length, file_percent) / rel)
                write_phylip(trim_sequences(reduced, length), dst)
                dsts.append(dst)
    return dsts


def run(input_dir: Path, output_dir: Path, lengths: Sequence[int],
        workers: int | None = None,
        missing: Sequence[float] = (),
        outgroup: str | None = DEFAULT_OUTGROUP,
        seed: int = DEFAULT_SEED,
        file_percents: Sequence[float] = ()) -> None:
    """Stream files from ``input_dir`` and process them in parallel.

    Files are never all held in memory: the directory walk is consumed
    lazily and only ``workers * 4`` tasks are in flight at any time. A quick
    metadata-only pre-walk counts the files up front so the progress bar can
    show a percentage.

    Progress reporting adapts to the output: a live tqdm bar when stderr is
    a terminal, and periodic one-line-per-update progress messages when the
    output is piped or redirected (where ``\\r``-based bars cannot render).
    """
    files = list(find_phy_files(input_dir))
    total = len(files)
    if total == 0:
        print(f"warning: no .phy files found under {input_dir}", file=sys.stderr)
        return

    if file_percents and not missing:
        raise ValueError(
            "file_percents requires missing_percents: with no taxa removed "
            "there is nothing to restrict to a subset of the files")

    # Threshold tables for the exact per-directory selection. Computed once
    # here rather than per worker, so the selected count is exact and the same
    # for any --workers.
    cutoffs = {f: selection_cutoffs(files, input_dir, seed, f)
               for f in file_percents if f < 100}

    if workers is None:
        workers = os.cpu_count() or 1
    workers = max(1, min(workers, total))

    levels = [f"{len(lengths)} lengths", f"{len(missing)} missing levels"]
    if file_percents and any(f != 100 for f in file_percents):
        levels.append(f"{len(file_percents)} file levels")
    desc = f"processing {total} files ({', '.join(levels)})"
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
            for src in files:
                process_file(src, input_dir, output_dir, lengths,
                             missing, outgroup, seed, file_percents, cutoffs)
                advance(1)
            return

        # Bounded window of in-flight futures over the already-listed files.
        max_in_flight = workers * 4
        file_iter: Iterator[Path] = iter(files)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            pending = {
                pool.submit(process_file, src, input_dir, output_dir, lengths,
                            missing, outgroup, seed, file_percents, cutoffs)
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
                                    lengths, missing, outgroup, seed,
                                    file_percents, cutoffs))
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
    parser.add_argument("--file-percents", type=float, nargs="+",
                        help="percentages of .phy files that get taxa removed; "
                             "crossed with --missing and --lengths "
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
        file_percents = [float(x) for x in cfg.get("file_percents", [])]
        outgroup = cfg.get("outgroup", DEFAULT_OUTGROUP)
        seed = cfg.get("seed", DEFAULT_SEED)
        workers = cfg.get("workers")
    elif args.input and args.output and args.lengths:
        input_dir, output_dir, lengths = args.input, args.output, args.lengths
        missing, file_percents = [], []
        outgroup, seed = DEFAULT_OUTGROUP, DEFAULT_SEED
        workers = None
    else:
        parser.error("provide either --config or all of --input/--output/--lengths")

    if args.missing is not None:
        missing = args.missing
    if args.file_percents is not None:
        file_percents = args.file_percents
    if args.outgroup is not None:
        outgroup = args.outgroup
    if args.seed is not None:
        seed = args.seed
    if args.workers is not None:
        workers = args.workers

    run(input_dir, output_dir, lengths, workers, missing, outgroup, seed,
        file_percents)


if __name__ == "__main__":
    main()
