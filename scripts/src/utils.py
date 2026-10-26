import re
import tomllib
import os
import pandas as pd

def remove_branch_lengths(net):
    x = re.sub(r":[0-9.:]*(?=[,)])", "", net)
    return re.sub(r"(?<=\))[0-9.:]+(?=[,;)])", "", x)


def _parse_enewick(net):
    """Parse an extended Newick network string into a node tree.

    Each node is a dict with keys:
      name       -- label (None if absent)
      hybrid     -- hybrid label without '#' (e.g. 'H31'), None if not hybrid
      is_ref     -- True for a hybrid reference leaf (bare '#Hxx')
      length     -- branch length of the edge into this node (None if absent)
      gamma      -- inheritance probability of the edge into this node
                     (None if absent)
      children   -- list of child nodes
    """
    net = net.strip().rstrip(";")
    pos = 0

    def parse_node():
        nonlocal pos
        children = []
        name = None
        hybrid = None
        is_ref = False
        length = None
        gamma = None
        if pos < len(net) and net[pos] == "(":
            pos += 1
            while True:
                children.append(parse_node())
                if pos >= len(net):
                    raise ValueError("Unbalanced parentheses in network")
                if net[pos] == ",":
                    pos += 1
                    continue
                if net[pos] == ")":
                    pos += 1
                    break
        match = re.match(r"#([A-Za-z0-9_]+)", net[pos:])
        if match:
            hybrid = match.group(1)
            pos += match.end()
            is_ref = not children
        match = re.match(r"[A-Za-z0-9_]+", net[pos:])
        if match:
            name = match.group(0)
            pos += match.end()
        if net[pos : pos + 1] == ":":
            pos += 1
            match = re.match(r"[-+0-9.eE]+", net[pos:])
            if match:
                length = float(match.group(0))
                pos += match.end()
            if net[pos : pos + 2] == "::":
                pos += 2
                match = re.match(r"[-+0-9.eE]+", net[pos:])
                if match:
                    gamma = float(match.group(0))
                    pos += match.end()
        return {
            "name": name,
            "hybrid": hybrid,
            "is_ref": is_ref,
            "length": length,
            "gamma": gamma,
            "children": children,
        }

    root = parse_node()
    if pos != len(net):
        raise ValueError(f"Trailing characters in network: {net[pos:]!r}")
    return root


def _collapse_unary(node, is_root=False):
    """Recursively collapse unary tree nodes (in-degree 1, out-degree 1).

    Mirrors PhyloNet's Networks.removeBinaryNodes: the two branch lengths
    are summed and the child's inheritance probability is preserved on the
    merged edge. Hybrid nodes (in-degree 2) and hybrid references are left
    alone; only plain tree nodes with exactly one child are collapsed.
    """
    node["children"] = [
        collapsed
        for child in node["children"]
        for collapsed in _collapse_unary(child)
    ]
    if not is_root and not node["hybrid"] and len(node["children"]) == 1:
        child = node["children"][0]
        # Sum the two branch lengths, treating a missing length as 1.0 --
        # the same default PhyloNet applies to the -s starting network before
        # its own (never-invoked here) binary-node removal would run.
        merged_length = (node["length"] if node["length"] is not None else 1.0) + (
            child["length"] if child["length"] is not None else 1.0
        )
        # The edge into the collapsed child keeps its own gamma (it may be a
        # hybrid edge); the gamma on the removed node's edge cannot exist
        # (a gamma requires the child to be hybrid, but this node is not).
        child["length"] = merged_length
        return [child]
    return [node]


def _format_number(value):
    if value is None:
        return None
    if value == int(value):
        return str(int(value))
    return repr(value)


def _to_enewick(node):
    """Render a node tree back to extended Newick."""
    parts = []
    if node["children"]:
        parts.append("(" + ",".join(_to_enewick(child) for child in node["children"]) + ")")
    if node["hybrid"]:
        parts.append("#" + node["hybrid"])
    elif node["name"]:
        parts.append(node["name"])
    length = _format_number(node["length"])
    if length is not None:
        parts.append(":" + length)
        if node["gamma"] is not None:
            parts.append("::" + _format_number(node["gamma"]))
    return "".join(parts)


def sanitize_network(net):
    """Sanitize a network for PhyloNet's InferNetwork_MPL (-s starting network).

    Fixes the structural defects that crash PhyloNet's search:
      * unary tree nodes (in-degree 1, out-degree 1), which make
        GeneTreeProbabilityPseudo throw NoSuchElementException and
        EdgeSourceChange throw 'should have two children'. These are
        collapsed, summing branch lengths and preserving hybrid gammas
        (same semantics as PhyloNet's own Networks.removeBinaryNodes).
      * a root with a single child is re-rooted on that child.

    Branch lengths and inheritance probabilities are preserved; PhyloNet
    keeps explicit values and only defaults missing ones (to 1.0 / 0.5).
    """
    root = _parse_enewick(net)
    root = _collapse_unary(root, is_root=True)[0]
    # Re-root if the root itself is unary (out-degree 1): the network is
    # equivalent to one rooted at its only child, and PhyloNet's search
    # assumes the root has two children.
    if len(root["children"]) == 1:
        only_child = root["children"][0]
        if only_child["length"] is not None:
            # The old root branch is dropped; the child keeps its own edge.
            pass
        root = only_child
        root["length"] = None
        root["gamma"] = None
    return _to_enewick(root) + ";"

def get_mpl_score(path):
    with open(path, "r") as file:
        for line_number, line in enumerate(file, start=1):
            if line_number == 5:
                return float(line.split(":", 1)[1].strip())

def get_fn_fp(path):
    with open(path, "r") as file:
        for line_number, line in enumerate(file, start=1):
            if line_number == 2:
                info = (line.split(":", 1)[1].strip()).split()
                return float(info[0]), float(info[1])

def get_fn_fp_fa(path):
    with open(path, "r") as file:
        for line_number, line in enumerate(file, start=1):
            if line_number == 2:
                info = (line.split(":", 1)[1].strip()).split()
                return float(info[0]), float(info[1]), float(info[2])

def load_config(path, required_keys):
    with open(path, "rb") as config_file:
        config = tomllib.load(config_file)

    missing_keys = required_keys - config.keys()

    if missing_keys:
        missing = ", ".join(sorted(missing_keys))
        raise ValueError(
            f"Configuration file is missing required keys: {missing}"
        )

    return config

def merge_records_with_stats(df_records, df_stats, how='inner'):
    """Merge accuracy metrics with performance metrics."""
    merged = df_records.merge(
        df_stats,
        left_on=['file_path', 'replicate'],
        right_on=['output_file', 'replicate'],
        how=how
    )
    # Drop the redundant output_file column
    merged = merged.drop('file_path', axis=1)
    return merged

def cmpnets_consolidate_records(input_path, file_names, method_labels):
    records = []

    for root, dirs, files in os.walk(input_path):
        replicate_name = os.path.relpath(root, input_path)
        for file_name in files:
            if file_name in file_names:
                file_path = os.path.join(root, file_name)
                fn, fp, fa = get_fn_fp_fa(file_path)
                records.append({
                    "method": file_name,
                    "file_path": file_path,
                    "false_negative": fn,
                    "false_positive": fp,
                    "false_average": fa,
                    "replicate": replicate_name,
                })

    df = pd.DataFrame(records)
    df["method_label"] = df["method"].map(method_labels)

    return df
