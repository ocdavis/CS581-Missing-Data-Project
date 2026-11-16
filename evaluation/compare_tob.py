#!/usr/bin/env python3
"""False negative / false positive branches of estimated trees vs. references.

Trees are compared as sets of non-trivial unrooted bipartitions (internal
branches), so rooting, branch lengths, support values and annotations are
ignored and polytomies are allowed:

  FN = branches of the reference missing from the estimate
  FP = branches of the estimate that are not in the reference
  fn_rate = FN / (internal branches of the reference)
  fp_rate = FP / (internal branches of the estimate)   (0 if it has none)

Two trees::

    python3 evaluation/compare_tob.py --estimated est.tre --reference ref.tre

Many trees: every file under ``--estimated`` whose name matches one of the
``--glob`` patterns is compared with its reference, looked up under
``--reference`` first in the same relative folder and then in a folder named
like the estimate's own folder (the replicate)::

    python3 evaluation/compare_tob.py --estimated data/gene_trees \\
        --glob tob-qmc.nwk "*.ectoblob.tre" \\
        --reference data/camus-dataset/n50 --reference-name true_tob.nwk \\
        --output results/errors.csv

One CSV row is written per estimated tree, and the mean rates per method
(file name) and condition (folders above the replicate) are printed.

By default the two trees must have the same taxa; ``--common-taxa`` instead
restricts both to the taxa they share. No third-party packages are needed.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

FIELDS = ["method", "condition", "replicate", "n_taxa", "ref_edges",
          "est_edges", "fn", "fp", "fn_rate", "fp_rate", "file", "reference"]


def leaf_sets(newick: str) -> tuple[frozenset[str], list[frozenset[str]]]:
    """Return (all leaves, leaf set below every internal node) of a tree.

    Handles quoted labels (which may contain any character), [comments],
    branch lengths, internal node labels and polytomies. Only the first tree
    of the text is read.
    """
    stack: list[list[str]] = []   # leaves collected for each open '('
    clades: list[frozenset[str]] = []
    leaves: list[str] = []
    after_close = False            # a label here names an internal node
    i, n = 0, len(newick)
    done = False
    while i < n and not done:
        c = newick[i]
        if c == "(":
            stack.append([])
            after_close = False
            i += 1
        elif c == ",":
            after_close = False
            i += 1
        elif c == ")":
            if not stack:
                raise ValueError("unbalanced Newick: unexpected ')'")
            below = stack.pop()
            clades.append(frozenset(below))
            if stack:
                stack[-1].extend(below)
            after_close = True
            i += 1
        elif c == ";":
            done = True
        elif c == "[":
            j = newick.find("]", i)
            if j < 0:
                raise ValueError("unterminated [comment]")
            i = j + 1
        elif c == ":":
            i += 1
            while i < n and newick[i] not in "(),;[":
                i += 1
        elif c.isspace():
            i += 1
        else:
            if c == "'":
                j, parts = i + 1, []
                while True:
                    k = newick.find("'", j)
                    if k < 0:
                        raise ValueError("unterminated quoted label")
                    parts.append(newick[j:k])
                    if newick[k + 1:k + 2] == "'":   # '' is an escaped quote
                        parts.append("'")
                        j = k + 2
                    else:
                        i = k + 1
                        break
                label = "".join(parts)
            else:
                j = i
                while j < n and newick[j] not in "(),:;[" and not newick[j].isspace():
                    j += 1
                label, i = newick[i:j], j
            if not after_close:
                if not stack:
                    raise ValueError("label outside parentheses")
                stack[-1].append(label)
                leaves.append(label)
    if stack or not done:
        raise ValueError("unbalanced or unterminated Newick string")
    if len(set(leaves)) != len(leaves):
        raise ValueError("duplicate taxon names")
    return frozenset(leaves), clades


def bipartitions(newick: str, keep: frozenset[str] | None = None
                 ) -> tuple[frozenset[str], set[frozenset[str]]]:
    """Return (taxa, non-trivial unrooted bipartitions), optionally after
    restricting the tree to the taxa in ``keep``.

    Each bipartition is stored as the side without the alphabetically
    smallest taxon, so sets from different trees are directly comparable.
    """
    taxa, clades = leaf_sets(newick)
    if keep is not None:
        taxa = taxa & keep
        clades = [c & keep for c in clades]
    if not taxa:
        return taxa, set()
    anchor = min(taxa)
    splits = set()
    for side in clades:
        if anchor in side:
            side = taxa - side
        if 2 <= len(side) <= len(taxa) - 2:
            splits.add(side)
    return taxa, splits


def compare(ref_newick: str, est_newick: str, common_taxa: bool = False) -> dict:
    ref_taxa, ref_splits = bipartitions(ref_newick)
    est_taxa, est_splits = bipartitions(est_newick)
    if ref_taxa != est_taxa:
        if not common_taxa:
            raise ValueError(
                f"different taxa ({len(ref_taxa - est_taxa)} only in the "
                f"reference, {len(est_taxa - ref_taxa)} only in the estimate); "
                f"use --common-taxa to compare on the shared ones")
        shared = ref_taxa & est_taxa
        if len(shared) < 4:
            raise ValueError("fewer than 4 taxa in common")
        ref_taxa, ref_splits = bipartitions(ref_newick, shared)
        est_taxa, est_splits = bipartitions(est_newick, shared)
    fn = len(ref_splits - est_splits)
    fp = len(est_splits - ref_splits)
    return {
        "n_taxa": len(ref_taxa),
        "ref_edges": len(ref_splits),
        "est_edges": len(est_splits),
        "fn": fn,
        "fp": fp,
        "fn_rate": fn / len(ref_splits) if ref_splits else 0.0,
        "fp_rate": fp / len(est_splits) if est_splits else 0.0,
    }


def find_reference(est: Path, est_root: Path, ref_root: Path, name: str) -> Path:
    """Reference for an estimate: same relative folder, else <replicate>/."""
    if ref_root.is_file():
        return ref_root
    mirrored = ref_root / est.parent.relative_to(est_root) / name
    if mirrored.exists():
        return mirrored
    return ref_root / est.parent.name / name


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        epilog="See the top of this file for examples.")
    parser.add_argument("--estimated", type=Path, required=True,
                        help="estimated tree file, or a directory to search")
    parser.add_argument("--reference", "--true", type=Path, required=True,
                        dest="reference",
                        help="reference tree file (used for every estimate), "
                             "or a directory holding one reference per folder")
    parser.add_argument("--glob", nargs="+", default=["tob-qmc.nwk"],
                        help="file name pattern(s) of the estimated trees when "
                             "--estimated is a directory (default: tob-qmc.nwk)")
    parser.add_argument("--reference-name", default="true_tob.nwk",
                        help="file name of the references when --reference is "
                             "a directory (default: true_tob.nwk)")
    parser.add_argument("--common-taxa", action="store_true",
                        help="restrict both trees to their shared taxa instead "
                             "of failing when the taxon sets differ")
    parser.add_argument("--output", type=Path, help="CSV file (default: stdout)")
    args = parser.parse_args()

    if args.estimated.is_file():
        est_root = args.estimated.parent
        files = [args.estimated]
    else:
        est_root = args.estimated
        files = sorted({p for g in args.glob for p in est_root.rglob(g) if p.is_file()})
        if not files:
            sys.exit(f"no files matching {args.glob} under {est_root}")

    rows, failed = [], 0
    for est in files:
        ref = find_reference(est, est_root, args.reference, args.reference_name)
        if ref.resolve() == est.resolve():
            continue  # a reference found by the same pattern
        try:
            row = compare(ref.read_text(), est.read_text(), args.common_taxa)
        except (OSError, ValueError) as exc:
            failed += 1
            print(f"ERROR: {est}: {exc}", file=sys.stderr)
            continue
        rel = est.relative_to(est_root)
        rows.append({
            "method": est.name,
            "condition": "/".join(rel.parts[:-2]),
            "replicate": rel.parts[-2] if len(rel.parts) > 1 else "",
            **row,
            "file": rel.as_posix(),
            "reference": ref.as_posix(),
        })

    out = args.output.open("w", newline="") if args.output else sys.stdout
    writer = csv.DictWriter(out, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({**row, "fn_rate": f"{row['fn_rate']:.4f}",
                         "fp_rate": f"{row['fp_rate']:.4f}"})
    if args.output:
        out.close()

    # Summary per method and condition (to stderr so stdout stays pure CSV).
    groups = defaultdict(list)
    for row in rows:
        groups[(row["method"], row["condition"])].append(row)
    log = sys.stderr
    print(f"\n{'method':<28}{'condition':<20}{'trees':>6}{'FN':>8}{'FP':>8}"
          f"{'FN rate':>9}{'FP rate':>9}", file=log)
    for (method, cond), g in sorted(groups.items()):
        k = len(g)
        print(f"{method:<28}{cond or '-':<20}{k:>6}"
              f"{sum(r['fn'] for r in g) / k:>8.2f}{sum(r['fp'] for r in g) / k:>8.2f}"
              f"{sum(r['fn_rate'] for r in g) / k:>9.4f}"
              f"{sum(r['fp_rate'] for r in g) / k:>9.4f}", file=log)
    if args.output:
        print(f"\nwrote {len(rows)} rows to {args.output}", file=log)
    if failed:
        sys.exit(f"{failed} of {len(files)} files could not be compared")


if __name__ == "__main__":
    main()
