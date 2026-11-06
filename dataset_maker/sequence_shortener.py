"""Create shortened copies of PHYLIP alignment files.

For every ``*.phy`` file under ``input_dir`` (recursively) and every length in the
config file, write a copy with only the first ``length`` characters of each
sequence.  Taxon labels and the file layout are kept as they are; only the
sequence characters and the length in the header line change.

Usage::

    python dataset_maker/seqence_shortener.py dataset_maker/shorten_example.toml

Config file (TOML)::

    input_dir  = "archives/camus-dataset/n25/03/seqs-500"
    output_dir = "data/shortened"
    lengths    = [100, 200, 300, 400]
    steps      = ["truncate"]                    # optional, this is the default

A file ``<input_dir>/n25/03/seqs-500/0002.phy`` then produces::

    <output_dir>/100/n25/03/seqs-500/0002.phy
    <output_dir>/200/n25/03/seqs-500/0002.phy
    ...

Adding a new operation
----------------------
Each file is read once into an :class:`Alignment` (a list of :class:`Row` objects
plus the sequence length).  The per-length work is a list of *steps* named in the
config and applied left to right.  To add one:

1. write a function ``(alignment, length, **options) -> Alignment``;
2. add it to :data:`STEPS`.

Example - drop some taxa, then in the config write
``steps = ["truncate", {action = "drop_taxa", taxa = ["t5", "t12"]}]``::

    def drop_taxa(alignment, length, *, taxa):
        keep = [row for row in alignment.rows if row.label not in set(taxa)]
        return Alignment(keep, alignment.nchar, alignment.tail)

    STEPS["drop_taxa"] = drop_taxa

A step must return a *new* ``Alignment`` instead of modifying the one it is
given, because the same parsed alignment is reused for every length.  The
``<ntax>`` in the header is always recomputed from the rows, so steps may add or
remove taxa freely.
"""

from __future__ import annotations

import re
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path

#: A data row is ``<name><whitespace><sequence>``; group 1 is the name column
#: (kept verbatim, padding included), group 2 is the sequence.
_ROW = re.compile(r"^(\S+\s+)(\S+)")

#: Name of the default step, so a config without ``steps`` still does something.
DEFAULT_STEPS = ("truncate",)


@dataclass
class Row:
    """One sequence line of a PHYLIP file."""

    column: str  # name column exactly as written, e.g. "OUT     "
    sequence: str  # sequence characters, e.g. "ACGTACGT..."

    @property
    def label(self) -> str:
        """Taxon name without its padding, e.g. ``"OUT"``."""
        return self.column.strip()


@dataclass
class Alignment:
    """A parsed file: its data rows, the sequence length, and the trailing text."""

    rows: tuple[Row, ...]
    nchar: int
    tail: str = "\n"  # whatever followed the last row (usually blank lines)

    def __post_init__(self) -> None:
        self.rows = tuple(self.rows)
        if not self.rows:
            raise ValueError("no sequences found")
        for row in self.rows:
            if len(row.sequence) != self.nchar:
                raise ValueError(
                    f"row {row.label!r} has {len(row.sequence)} sites, "
                    f"expected {self.nchar}"
                )

    @property
    def ntax(self) -> int:
        """Number of rows (recomputed, so steps may add or remove taxa)."""
        return len(self.rows)


def read_alignment(path: Path) -> Alignment:
    """Parse a sequential PHYLIP file (one sequence per line).

    Blank lines, trailing blank lines and the exact name column are preserved, so
    reading then writing a file back changes nothing except the padding at the end
    of a data row (the CAMUS files pad rows with spaces; those are dropped because
    no PHYLIP reader looks at them).
    """
    lines = path.read_text().split("\n")

    header_index = next((i for i, line in enumerate(lines) if line.strip()), None)
    if header_index is None:
        raise ValueError("file is empty")

    fields = lines[header_index].split()
    if len(fields) < 2 or not (fields[0].isdigit() and fields[1].isdigit()):
        raise ValueError(f"bad header line: {lines[header_index]!r}")

    rows: list[Row] = []
    last_row = header_index
    for index in range(header_index + 1, len(lines)):
        match = _ROW.match(lines[index])
        if match:
            rows.append(Row(match.group(1), match.group(2)))
            last_row = index
        elif lines[index].strip():
            raise ValueError(
                f"line {index + 1}: expected '<name> <sequence>', "
                f"got {lines[index]!r}"
            )

    tail_lines = lines[last_row + 1 :]
    tail = "\n" + "\n".join(tail_lines) if tail_lines else "\n"
    return Alignment(rows, int(fields[1]), tail)


def format_alignment(alignment: Alignment) -> str:
    """Render an :class:`Alignment` back to PHYLIP text."""
    lines = [f"{alignment.ntax}  {alignment.nchar}"]
    lines.extend(f"{row.column}{row.sequence}" for row in alignment.rows)
    return "\n".join(lines) + alignment.tail


def write_alignment(alignment: Alignment, target: Path) -> None:
    """Write ``alignment`` to ``target``, creating parent directories as needed."""
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(format_alignment(alignment))


# --- steps -------------------------------------------------------------------


def truncate(alignment: Alignment, length: int, **_options) -> Alignment:
    """Keep the first ``length`` sites of every sequence."""
    if length > alignment.nchar:
        raise ValueError(f"alignment is only {alignment.nchar} sites long")
    rows = [Row(row.column, row.sequence[:length]) for row in alignment.rows]
    return Alignment(rows, length, alignment.tail)


#: Registry of steps, keyed by the ``action`` names used in the config file.
STEPS = {"truncate": truncate}


def apply_steps(alignment: Alignment, steps, length: int) -> Alignment:
    """Apply the configured steps, left to right, for one target length.

    A step is either a name (``"truncate"``) or a table with an ``action`` key
    plus options (``{action = "drop_taxa", taxa = ["t5"]}``).
    """
    result = alignment
    for step in steps:
        if isinstance(step, dict):
            options = dict(step)
            action = options.pop("action", None)
        else:
            action, options = step, {}

        try:
            function = STEPS[action]
        except KeyError:
            raise ValueError(
                f"unknown step {action!r}; available: {sorted(STEPS)}"
            ) from None
        result = function(result, length, **options)
    return result


# --- config and driver -------------------------------------------------------


def load_config(path: Path) -> dict:
    """Read the TOML config; relative paths are resolved against the config file."""
    with path.open("rb") as handle:
        raw = tomllib.load(handle)

    for key in ("input_dir", "output_dir", "lengths"):
        if key not in raw:
            raise ValueError(f"{path}: missing {key!r}")
    if not isinstance(raw["lengths"], list) or not raw["lengths"]:
        raise ValueError(f"{path}: 'lengths' must be a non-empty list")

    steps = raw.get("steps", list(DEFAULT_STEPS))
    if not isinstance(steps, list) or not steps:
        raise ValueError(f"{path}: 'steps' must be a non-empty list")

    base = path.resolve().parent
    return {
        "input_dir": (base / raw["input_dir"]).resolve(),
        "output_dir": (base / raw["output_dir"]).resolve(),
        "lengths": sorted({int(length) for length in raw["lengths"]}),
        "steps": steps,
    }


def build_datasets(
    input_dir: Path,
    output_dir: Path,
    lengths: list[int],
    steps=DEFAULT_STEPS,
) -> None:
    """Write one dataset per length, mirroring the input directory structure."""
    sources = sorted(path for path in input_dir.rglob("*.phy") if path.is_file())
    print(f"{len(sources)} file(s) found in {input_dir}")

    written = 0
    for source in sources:
        relative = source.relative_to(input_dir)
        try:
            alignment = read_alignment(source)
        except (OSError, ValueError) as exc:
            print(f"skipping {source}: {exc}", file=sys.stderr)
            continue

        for length in lengths:
            try:
                result = apply_steps(alignment, steps, length)
            except (TypeError, ValueError) as exc:
                print(f"skipping {source} ({length}): {exc}", file=sys.stderr)
                continue
            write_alignment(result, output_dir / str(length) / relative)
            written += 1

    print(f"{written} file(s) written to {output_dir}")


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(f"usage: {Path(sys.argv[0]).name} CONFIG.toml", file=sys.stderr)
        return 2

    try:
        config = load_config(Path(argv[0]))
    except (OSError, ValueError, tomllib.TOMLDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    build_datasets(**config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
