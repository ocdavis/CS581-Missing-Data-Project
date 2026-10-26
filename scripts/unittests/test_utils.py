import unittest
from src.utils import *
from src.utils import _parse_enewick
import sys
import os
sys.dont_write_bytecode = True


class TestUtils(unittest.TestCase):

    def test_remove_branch_lengths(self):
        input1 = "((t1:0.12029,t2:0.232309)(t3:23.237323,t4:893.829389));"
        expected_output1 = "((t1,t2)(t3,t4));"
        calculated_output1 = remove_branch_lengths(input1)
        self.assertEqual(expected_output1, calculated_output1)
        input2 = "(t3,t12,(t2,((((t13,(t11,t7)1:0.33894013053783395)1:0.45718108964384196,((t15,t10):0.021230360592972005,#H30:0.04634541896266962::0.3539320314293897):1.0721775005863046)1:0.48235201057096455,((t1,t9)1:0.040234728770749704,((t14,t4)0.985724:0.054416203199337806,(t5,((t6,t8)1:0.0)#H30:0.08039797772894099::0.6460679685706103)0.532078:0.012528005514031309)1:0.3308743841186229)0.997299:1.3089997541445408)1:0.4885065344542067,OUT)1:1.3011012767550771)1:0.8502460402338993)1;"
        expected_output2 = "(t3,t12,(t2,((((t13,(t11,t7)),((t15,t10),#H30)),((t1,t9),((t14,t4),(t5,((t6,t8))#H30)))),OUT)));"
        calculated_output2 = remove_branch_lengths(input2)
        self.assertEqual(expected_output2, calculated_output2)

    def test_get_mpl_score(self):
        loc_path = os.path.dirname(__file__)
        expected_mpl = -487565.7058449446
        calculated_mpl = get_mpl_score(os.path.join(loc_path, "calgtprob_example.txt"))
        self.assertEqual(expected_mpl, calculated_mpl)

    def test_get_fn_fp(self):
        loc_path = os.path.dirname(__file__)
        expected_fn = 0.21052631578947367
        expected_fp = 0.11764705882352941
        calculated_fn, calculated_fp = get_fn_fp(os.path.join(loc_path, "cmpnets_example.txt"))
        self.assertEqual(expected_fn, calculated_fn)
        self.assertEqual(expected_fp, calculated_fp)


def _traverse(node, is_root=True):
    """Recursively collect (node_key, children_keys) pairs for topology comparison.

    Hybrid nodes are keyed by their '#Hxx' label so that a hybrid definition
    and its reference are recognized as the same node; plain tree nodes are
    keyed by their sorted child keys (unary nodes thus vanish naturally).
    """
    if node["is_ref"]:
        return "#" + node["hybrid"]
    child_keys = [_traverse(c, is_root=False) for c in node["children"]]
    if node["hybrid"]:
        return "#" + node["hybrid"] + "(" + ",".join(child_keys) + ")"
    if node["name"] and not child_keys:
        return node["name"]
    if len(child_keys) == 1:
        # Unary wrapper (including the root): transparent for topology --
        # pass through its child.
        return child_keys[0]
    return "(" + ",".join(child_keys) + ")"


def _topology_key(net):
    """Canonical topology signature: ignores branch lengths, gammas, and
    unary wrappers, but keeps leaf labels, cluster structure, and the
    reticulation structure (hybrid labels and their placements)."""
    return _traverse(_parse_enewick(net))


def _hybrid_definitions(node, defs):
    """Map each hybrid label to its definition node."""
    if node["hybrid"] and not node["is_ref"]:
        defs[node["hybrid"]] = node
    for c in node["children"]:
        _hybrid_definitions(c, defs)
    return defs


def _leaf_set(node, acc, defs=None, _seen=None):
    """Collect leaf labels, resolving hybrid references to the leaves under
    the hybrid's definition node."""
    if defs is None:
        defs = _hybrid_definitions(node, {})
        _seen = set()
    if node["is_ref"]:
        if node["hybrid"] in _seen:  # avoid infinite recursion on cycles
            return acc
        _seen.add(node["hybrid"])
        return _leaf_set(defs[node["hybrid"]], acc, defs, _seen)
    if node["name"] and not node["children"]:
        acc.add(node["name"])
    for c in node["children"]:
        _leaf_set(c, acc, defs, _seen)
    return acc


def _clusters(node, acc, defs=None, _seen=None):
    """Collect all clusters (resolved leaf sets under each tree node),
    treating hybrid references as their definition's leaf set."""
    if node["is_ref"]:
        return acc
    if defs is None:
        defs = _hybrid_definitions(node, {})
        _seen = set()
    if node["children"] and not node["hybrid"]:
        leaves = set()
        for c in node["children"]:
            _leaf_set(c, leaves, defs, set(_seen))
        acc.append(frozenset(leaves))
    for c in node["children"]:
        _clusters(c, acc, defs, _seen)
    return acc


class TestSanitizeNetwork(unittest.TestCase):

    def test_unary_node_under_root_is_collapsed(self):
        """CAMUS root artifact: (OUT,(BIG)len)len -- the (BIG) wrapper must go."""
        net = "(OUT:0.9,(((t1:0.5,t2:0.5):0.6,(t3:0.7,t4:0.7):0.8):0.2):0.1);"
        sanitized = sanitize_network(net)
        self.assertEqual(_topology_key(sanitized), _topology_key(net))
        # The unary wrapper is gone: root's second child is the subtree itself.
        self.assertEqual(
            sanitized,
            "(OUT:0.9,((t1:0.5,t2:0.5):0.6,(t3:0.7,t4:0.7):0.8):0.30000000000000004);",
        )

    def test_topology_unchanged_on_hybrid_network(self):
        """Full hybrid network: topology must be identical after sanitizing."""
        net = "(OUT:0.9844421851525048,((((((t7:0.4268367413,t11:0.4268367413):0.3032192815,t13:0.7300560228):0.5570698435,(t10:0.1042706784,(t15:0.09129985249,#H31:0.0::0.3543951263):0.01297082592):1.182855188):0.5360690747,(((t5:0.09129985249,((t6:0.06460095234,t8:0.06460095234):0.02669890015)#H31:0.0::0.6456048737):0.02223104008,(t4:0.04046950944,t14:0.04046950944):0.07306138313):0.3255581916,(t1:0.3891996428,t9:0.3891996428):0.04988944144):1.384105857):0.6086537933,(t2:0.9483343827,(t3:0.01998138516,t12:0.01998138516):0.9283529975):1.483514352):1.33092468):0.07579544029403025);"
        sanitized = sanitize_network(net)
        # Same topology: identical canonical signature.
        self.assertEqual(_topology_key(sanitized), _topology_key(net))
        # Same clusters: every split of the original is present.
        self.assertEqual(sorted(_clusters(_parse_enewick(sanitized), [])),
                         sorted(_clusters(_parse_enewick(net), [])))
        # Same leaves.
        self.assertEqual(_leaf_set(_parse_enewick(sanitized), set()),
                         _leaf_set(_parse_enewick(net), set()))
        # Hybrid edges and gammas preserved.
        self.assertIn("#H31:0::0.3543951263", sanitized)
        self.assertIn("#H31:0::0.6456048737", sanitized)

    def test_no_unary_nodes_remain(self):
        net = "(OUT:0.9,(((t1:0.5,t2:0.5):0.6,(t3:0.7,t4:0.7):0.8):0.2):0.1);"
        sanitized = sanitize_network(net)

        def count_unary(node, is_root=True):
            n = 0
            if (not is_root and not node["is_ref"] and not node["hybrid"]
                    and len(node["children"]) == 1):
                n += 1
            for c in node["children"]:
                n += count_unary(c, is_root=False)
            return n

        self.assertEqual(count_unary(_parse_enewick(sanitized)), 0)

    def test_root_with_single_child_is_rerooted(self):
        """A root with out-degree 1 is re-rooted on its only child."""
        net = "(((t1:0.5,t2:0.5):0.6,(t3:0.7,t4:0.7):0.8):0.2);"
        sanitized = sanitize_network(net)
        self.assertEqual(_topology_key(sanitized), _topology_key(net))
        root = _parse_enewick(sanitized)
        self.assertEqual(len(root["children"]), 2)

    def test_already_clean_network_is_unchanged(self):
        net = "(OUT:0.9,((t1:0.5,t2:0.5):0.6,(t3:0.7,t4:0.7):0.8));"
        sanitized = sanitize_network(net)
        self.assertEqual(_topology_key(sanitized), _topology_key(net))
        # No unary wrappers existed, so the string should be identical.
        self.assertEqual(sanitized, net)

    def test_missing_lengths_defaulted_on_collapse(self):
        """A lengthless edge collapsed with a lengthed one sums as 1.0 + len."""
        net = "(OUT,(((t1,t2):0.6,(t3,t4):0.8):0.2));"
        sanitized = sanitize_network(net)
        # 1.0 (default for the missing outer edge) + 0.2 = 1.2 on the merged edge.
        self.assertIn(":1.2", sanitized)

    def test_malformed_network_raises(self):
        with self.assertRaises(ValueError):
            sanitize_network("(t1,t2")  # unbalanced

    def test_real_n15_networks_preserve_topology(self):
        """Every rooted n15 network: sanitized topology == original topology."""
        import glob
        base = os.path.join(os.path.dirname(__file__), "..", "..", "data",
                            "rooted-camus-dataset", "n15")
        files = sorted(glob.glob(os.path.join(base, "*", "true_net-rooted.nwk")))
        self.assertGreater(len(files), 0, "no n15 networks found")
        for f in files:
            with open(f) as fh:
                net = fh.read().strip()
            sanitized = sanitize_network(net)
            self.assertEqual(
                _topology_key(sanitized), _topology_key(net),
                msg=f"topology changed for {f}")
            self.assertEqual(
                sorted(_clusters(_parse_enewick(sanitized), [])),
                sorted(_clusters(_parse_enewick(net), [])),
                msg=f"clusters changed for {f}")
            self.assertEqual(
                _leaf_set(_parse_enewick(sanitized), set()),
                _leaf_set(_parse_enewick(net), set()),
                msg=f"leaf set changed for {f}")

