#!/usr/bin/env python3
"""Polytomies of estimated trees of blobs compared with the true ones.

A polytomy is an internal node of degree 4 or more in the unrooted tree; in
a tree of blobs each polytomy stands for a blob, and its degree is the number
of subtrees attached to that blob. For every estimated tree and its reference
this reports the number of polytomies and the largest polytomy degree (3 if
the tree is fully resolved).

    python3 evaluation/polytomy_stats.py --estimated data/gene_trees \\
        --glob tob-qmc.nwk "*.ectoblob.tre" \\
        --reference data/camus-dataset/n50 --output results/polytomy_stats.csv

Files are found exactly as in compare_tob.py (same options). One CSV row is
written per estimated tree and the means per method and condition are
printed. No third-party packages are needed.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from compare_tob import find_reference

FIELDS = ["method", "condition", "replicate", "n_taxa",
          "est_polytomies", "ref_polytomies", "est_max_degree",
          "ref_max_degree", "est_internal_nodes", "ref_internal_nodes",
          "file", "reference"]


def node_degrees(newick: str) -> tuple[int, list[int]]:
    """Return (number of leaves, unrooted degree of every internal node).

    A root with two children is not a node of the unrooted tree and is
    skipped. Quoted labels, [comments], branch lengths and internal labels
    are ignored; only the first tree of the text is read.
    """
    text = re.sub(r"'(?:[^']|'')*'", "", newick)      # quoted labels
    text = re.sub(r"\[[^\]]*\]", "", text)            # comments
    text = text.split(";")[0]
    stack: list[int] = []      # number of children seen for each open '('
    degrees: list[int] = []
    leaves = 0
    prev = ""
    for token in re.findall(r"[(),]|[^(),]+", text):
        if token == "(":
            if stack:
                stack[-1] += 1
            stack.append(0)
        elif token == ")":
            if not stack:
                raise ValueError("unbalanced Newick: unexpected ')'")
            children = stack.pop()
            if stack:
                degrees.append(children + 1)          # + the edge to its parent
            elif children != 2:
                degrees.append(children)              # a real root node
        elif token != "," and prev in ("(", ","):
            if token.split(":")[0].strip():           # a leaf label
                if not stack:
                    raise ValueError("label outside parentheses")
                stack[-1] += 1
                leaves += 1
        prev = token
    if stack:
        raise ValueError("unbalanced Newick string")
    return leaves, degrees


def stats(newick: str) -> dict:
    leaves, degrees = node_degrees(newick)
    return {"n_taxa": leaves,
            "polytomies": sum(d >= 4 for d in degrees),
            "max_degree": max(degrees, default=0),
            "internal_nodes": len(degrees)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--estimated", type=Path, required=True,
                        help="estimated tree file, or a directory to search")
    parser.add_argument("--reference", "--true", type=Path, required=True,
                        dest="reference",
                        help="reference tree file, or a directory holding "
                             "one reference per folder")
    parser.add_argument("--glob", nargs="+", default=["tob-qmc.nwk"],
                        help="file name pattern(s) of the estimated trees "
                             "(default: tob-qmc.nwk)")
    parser.add_argument("--reference-name", default="true_tob.nwk",
                        help="file name of the references (default: "
                             "true_tob.nwk)")
    parser.add_argument("--output", type=Path, help="CSV file (default: stdout)")
    args = parser.parse_args()

    if args.estimated.is_file():
        est_root, files = args.estimated.parent, [args.estimated]
    else:
        est_root = args.estimated
        files = sorted({p for g in args.glob for p in est_root.rglob(g)
                        if p.is_file()})
        if not files:
            sys.exit(f"no files matching {args.glob} under {est_root}")

    rows, failed = [], 0
    for est in files:
        ref = find_reference(est, est_root, args.reference, args.reference_name)
        try:
            e, r = stats(est.read_text()), stats(ref.read_text())
        except (OSError, ValueError) as exc:
            failed += 1
            print(f"ERROR: {est}: {exc}", file=sys.stderr)
            continue
        rel = est.relative_to(est_root)
        rows.append({
            "method": est.name,
            "condition": "/".join(rel.parts[:-2]),
            "replicate": rel.parts[-2] if len(rel.parts) > 1 else "",
            "n_taxa": e["n_taxa"],
            "est_polytomies": e["polytomies"],
            "ref_polytomies": r["polytomies"],
            "est_max_degree": e["max_degree"],
            "ref_max_degree": r["max_degree"],
            "est_internal_nodes": e["internal_nodes"],
            "ref_internal_nodes": r["internal_nodes"],
            "file": rel.as_posix(),
            "reference": ref.as_posix(),
        })

    out = args.output.open("w", newline="") if args.output else sys.stdout
    writer = csv.DictWriter(out, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    if args.output:
        out.close()

    groups = defaultdict(list)
    for row in rows:
        groups[(row["method"], row["condition"])].append(row)
    log = sys.stderr
    print(f"\n{'method':<26}{'condition':<28}{'trees':>6}"
          f"{'polytomies':>12}{'(true)':>8}{'max degree':>12}{'(true)':>8}",
          file=log)
    for (method, cond), g in sorted(groups.items()):
        k = len(g)
        mean = lambda col: sum(r[col] for r in g) / k   # noqa: E731
        print(f"{method:<26}{cond or '-':<28}{k:>6}"
              f"{mean('est_polytomies'):>12.2f}{mean('ref_polytomies'):>8.2f}"
              f"{mean('est_max_degree'):>12.2f}{mean('ref_max_degree'):>8.2f}",
              file=log)
    if args.output:
        print(f"\nwrote {len(rows)} rows to {args.output}", file=log)
    if failed:
        sys.exit(f"{failed} of {len(files)} files could not be read")


if __name__ == "__main__":
    main()
