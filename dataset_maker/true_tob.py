#!/usr/bin/env python3
"""Compute the true tree of blobs (ToB) of phylogenetic networks.

The tree of blobs of a network is obtained by collapsing every blob (a
maximal set of nodes held together by reticulation cycles) into a single
node. On the undirected graph of the network this is:

  1. find the cut edges (bridges) -- edges whose removal disconnects the
     graph; every edge that lies on a cycle is a non-bridge;
  2. contract all non-bridge edges, so each blob becomes one node;
  3. suppress nodes left with degree 2 (summing their branch lengths).

The result is a tree on the same leaves whose polytomies are the blobs.
Branch lengths of the surviving cut edges are kept. The output is written
rooted where the network was rooted, but for comparisons it should be
treated as unrooted.

For every network file found under the input directory a ToB file is
written next to it (or under ``--output``, mirroring the directory layout)::

    python3 true_tob.py --input data/camus-dataset/n50
    python3 true_tob.py --input data/camus-dataset/n50 --output data/true_tob

With ``--reference DIR`` each computed ToB is also compared, by its set of
unrooted bipartitions, against the ToB file at the same relative path under
DIR (e.g. the ``true_tob.nwk`` files distributed with BROOQS-Data).

No third-party packages are needed. ``tree_of_blobs`` and ``bipartitions``
can be imported and reused, e.g. for FN/FP error rates.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.setrecursionlimit(100_000)

# ---------------------------------------------------------------------------
# (Extended) Newick parsing
# ---------------------------------------------------------------------------


@dataclass
class Node:
    label: str = ""
    length: str | None = None  # kept as text so values round-trip exactly
    parent: "Node | None" = None
    children: list["Node"] = field(default_factory=list)


def parse_newick(text: str) -> Node:
    """Parse a Newick / extended Newick string into a tree of ``Node``s.

    Hybrid nodes (labels starting with ``#``) appear more than once in
    extended Newick; they are returned as separate ``Node``s here and merged
    later by ``tree_of_blobs``. Of ``label:length:support:gamma`` only the
    label and the length are kept.
    """
    text = text.strip()
    if not text.endswith(";"):
        raise ValueError("Newick string does not end with ';'")
    root = cur = Node()
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c == "(":
            child = Node(parent=cur)
            cur.children.append(child)
            cur = child
            i += 1
        elif c == ",":
            if cur.parent is None:
                raise ValueError("unbalanced Newick: ',' outside parentheses")
            sibling = Node(parent=cur.parent)
            cur.parent.children.append(sibling)
            cur = sibling
            i += 1
        elif c == ")":
            if cur.parent is None:
                raise ValueError("unbalanced Newick: unexpected ')'")
            cur = cur.parent
            i += 1
        elif c == ";":
            break
        elif c.isspace():
            i += 1
        else:
            j = i
            while j < n and text[j] not in "(),;":
                j += 1
            label, _, rest = text[i:j].partition(":")
            cur.label = label.strip()
            cur.length = rest.split(":")[0].strip() or None
            i = j
    if cur is not root:
        raise ValueError("unbalanced Newick: missing ')'")
    return root


def _walk(root: Node):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(node.children)


# ---------------------------------------------------------------------------
# Tree of blobs
# ---------------------------------------------------------------------------


def _find_bridges(n_nodes: int, edges: list[tuple[int, int]]) -> set[int]:
    """Return the indices of the bridge edges of an undirected multigraph."""
    adj: list[list[tuple[int, int]]] = [[] for _ in range(n_nodes)]
    for eid, (u, v) in enumerate(edges):
        adj[u].append((v, eid))
        adj[v].append((u, eid))
    disc = [-1] * n_nodes
    low = [0] * n_nodes
    bridges: set[int] = set()
    timer = 0
    for start in range(n_nodes):
        if disc[start] != -1:
            continue
        disc[start] = low[start] = timer
        timer += 1
        stack = [(start, -1, iter(adj[start]))]  # (node, edge we came by, it)
        while stack:
            u, via, it = stack[-1]
            for v, eid in it:
                if eid == via:  # skip only this edge, so parallel edges count
                    continue
                if disc[v] == -1:
                    disc[v] = low[v] = timer
                    timer += 1
                    stack.append((v, eid, iter(adj[v])))
                    break
                low[u] = min(low[u], disc[v])
            else:
                stack.pop()
                if stack:
                    parent = stack[-1][0]
                    low[parent] = min(low[parent], low[u])
                    if low[u] > disc[parent]:
                        bridges.add(via)
    return bridges


def _add_lengths(a: str | None, b: str | None) -> str | None:
    if a is None or b is None:
        return a if b is None else b
    return f"{float(a) + float(b):.10g}"


def tree_of_blobs(network_newick: str) -> str:
    """Return the tree of blobs of an extended Newick network, as Newick."""
    root = parse_newick(network_newick)

    # 1. Undirected graph; all occurrences of a hybrid label share one id.
    ids: dict[int, int] = {}
    hybrid_ids: dict[str, int] = {}
    leaf_name: dict[int, str] = {}
    n_nodes = 0
    nodes = list(_walk(root))
    for node in nodes:
        if node.label.startswith("#"):
            if node.label not in hybrid_ids:
                hybrid_ids[node.label] = n_nodes
                n_nodes += 1
            ids[id(node)] = hybrid_ids[node.label]
        else:
            ids[id(node)] = n_nodes
            if not node.children:
                if not node.label:
                    raise ValueError("unlabelled leaf in network")
                leaf_name[n_nodes] = node.label
            n_nodes += 1
    edges = [(ids[id(node.parent)], ids[id(node)]) for node in nodes
             if node.parent is not None]
    lengths = [node.length for node in nodes if node.parent is not None]

    # 2. Contract every non-bridge edge (union-find over blobs).
    bridges = _find_bridges(n_nodes, edges)
    comp = list(range(n_nodes))

    def find(x: int) -> int:
        while comp[x] != x:
            comp[x] = comp[comp[x]]
            x = comp[x]
        return x

    for eid, (u, v) in enumerate(edges):
        if eid not in bridges:
            comp[find(u)] = find(v)

    adj: dict[int, list[tuple[int, int]]] = {}
    for eid in bridges:
        u, v = find(edges[eid][0]), find(edges[eid][1])
        adj.setdefault(u, []).append((v, eid))
        adj.setdefault(v, []).append((u, eid))
    blob_leaf = {find(node): name for node, name in leaf_name.items()}

    # 3. Write the tree from the root's blob, suppressing degree-2 nodes.
    def build(blob: int, via: int) -> tuple[str, str | None]:
        """Newick of the subtree at ``blob`` and the length of the edge above."""
        length = lengths[via] if via >= 0 else None
        below = sorted((eid, nbr) for nbr, eid in adj.get(blob, []) if eid != via)
        if not below:
            if blob not in blob_leaf:
                raise ValueError("internal node without descendants in network")
            return blob_leaf[blob], length
        parts = [build(nbr, eid) for eid, nbr in below]
        if len(parts) == 1:
            sub, sub_length = parts[0]
            return sub, _add_lengths(sub_length, length)
        inner = ",".join(s if l is None else f"{s}:{l}" for s, l in parts)
        return f"({inner})", length

    newick, _ = build(find(ids[id(root)]), -1)
    return newick + ";"


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------


def bipartitions(tree_newick: str) -> tuple[frozenset[str], set[frozenset[str]]]:
    """Return ``(leaf set, non-trivial unrooted bipartitions)`` of a tree.

    Each bipartition is represented by the side that does not contain the
    alphabetically smallest leaf, so sets from two trees on the same leaves
    can be compared directly regardless of rooting.
    """
    root = parse_newick(tree_newick)
    below: dict[int, frozenset[str]] = {}

    def collect(node: Node) -> frozenset[str]:
        if not node.children:
            result = frozenset([node.label])
        else:
            result = frozenset().union(*(collect(c) for c in node.children))
        below[id(node)] = result
        return result

    leaves = collect(root)
    anchor = min(leaves)
    splits = set()
    for node in _walk(root):
        side = below[id(node)]
        if anchor in side:
            side = leaves - side
        if 2 <= len(side) <= len(leaves) - 2:
            splits.add(side)
    return leaves, splits


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", type=Path, required=True,
                        help="directory searched recursively for network files")
    parser.add_argument("--output", type=Path,
                        help="output directory mirroring the input layout "
                             "(default: write next to each network file)")
    parser.add_argument("--net-name", default="true_net.nwk",
                        help="network file name to look for (default: %(default)s)")
    parser.add_argument("--tob-name", default="true_tob.nwk",
                        help="ToB file name to write (default: %(default)s)")
    parser.add_argument("--reference", type=Path,
                        help="directory with reference ToB files (same relative "
                             "paths) to compare against")
    args = parser.parse_args(argv)

    nets = sorted(args.input.rglob(args.net_name))
    if not nets:
        print(f"error: no {args.net_name} files under {args.input}",
              file=sys.stderr)
        return 1

    out_root = args.output or args.input
    mismatches = compared = 0
    for net in nets:
        rel = net.parent.relative_to(args.input)
        tob = tree_of_blobs(net.read_text())
        dst = out_root / rel / args.tob_name
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(tob + "\n")

        leaves, splits = bipartitions(tob)
        line = (f"{rel.as_posix() or '.'}: {len(leaves)} leaves, "
                f"{len(splits)} internal branches "
                f"({len(leaves) - 3 - len(splits)} fewer than binary)")
        if args.reference:
            ref = args.reference / rel / args.tob_name
            if not ref.exists():
                line += " | no reference file"
            else:
                compared += 1
                ref_leaves, ref_splits = bipartitions(ref.read_text())
                if ref_leaves == leaves and ref_splits == splits:
                    line += " | matches reference"
                else:
                    mismatches += 1
                    line += (f" | DIFFERS from reference "
                             f"({len(splits - ref_splits)} extra, "
                             f"{len(ref_splits - splits)} missing)")
        print(line)

    print(f"wrote {len(nets)} ToB files under {out_root}")
    if args.reference:
        print(f"reference check: {compared - mismatches}/{compared} match")
    return 1 if mismatches else 0


if __name__ == "__main__":
    sys.exit(main())
