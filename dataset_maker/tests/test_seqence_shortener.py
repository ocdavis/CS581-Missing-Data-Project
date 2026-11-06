"""Tests for seqence_shortener.py.

Run with::

    python -m unittest discover -s dataset_maker/tests -t .
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from dataset_maker.sequence_shortener import (
    STEPS,
    Alignment,
    Row,
    apply_steps,
    build_datasets,
    format_alignment,
    load_config,
    read_alignment,
    truncate,
    write_alignment,
)

SAMPLE = """\
3  12
OUT     ACGTACGTACGT
t1      TTTTCCCCAAAA
t2      GGGGAAAATTTT

"""


class TempDirTestCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def sample_file(self, name: str = "0002.phy") -> Path:
        path = self.tmp / name
        path.write_text(SAMPLE)
        return path


class ReadAndWriteTests(TempDirTestCase):
    def test_round_trip_is_byte_identical(self):
        path = self.sample_file()
        self.assertEqual(format_alignment(read_alignment(path)), SAMPLE)

    def test_parses_rows_and_header(self):
        alignment = read_alignment(self.sample_file())
        self.assertEqual(alignment.ntax, 3)
        self.assertEqual(alignment.nchar, 12)
        self.assertEqual([row.label for row in alignment.rows], ["OUT", "t1", "t2"])
        self.assertEqual(alignment.rows[1].sequence, "TTTTCCCCAAAA")

    def test_rejects_a_file_without_a_phylip_header(self):
        path = self.tmp / "bad.phy"
        path.write_text("not a phylip file\n")
        with self.assertRaises(ValueError):
            read_alignment(path)

    def test_rejects_rows_of_different_lengths(self):
        path = self.tmp / "bad.phy"
        path.write_text("2  4\nOUT     ACGT\nt1      ACG\n")
        with self.assertRaisesRegex(ValueError, "has 3 sites, expected 4"):
            read_alignment(path)

    def test_write_creates_parent_directories(self):
        alignment = read_alignment(self.sample_file())
        target = self.tmp / "deep" / "nested" / "0002.phy"
        write_alignment(alignment, target)
        self.assertEqual(target.read_text(), SAMPLE)


class TruncateStepTests(TempDirTestCase):
    def test_keeps_first_sites_labels_and_layout(self):
        result = truncate(read_alignment(self.sample_file()), 5)
        lines = format_alignment(result).split("\n")
        self.assertEqual(lines[0], "3  5")
        self.assertEqual(lines[1], "OUT     ACGTA")
        self.assertEqual(lines[2], "t1      TTTTC")
        self.assertEqual(lines[3], "t2      GGGGA")
        self.assertTrue(format_alignment(result).endswith("t2      GGGGA\n\n"))

    def test_does_not_modify_the_input_alignment(self):
        original = read_alignment(self.sample_file())
        truncate(original, 5)
        self.assertEqual(original.nchar, 12)
        self.assertEqual(original.rows[0].sequence, "ACGTACGTACGT")

    def test_rejects_a_length_longer_than_the_alignment(self):
        alignment = read_alignment(self.sample_file())
        with self.assertRaisesRegex(ValueError, "only 12 sites"):
            truncate(alignment, 13)


class StepRegistryTests(TempDirTestCase):
    """The extension point: a human adds a function and one registry line."""

    def setUp(self):
        super().setUp()

        def drop_taxa(alignment, length, *, taxa):
            keep = [row for row in alignment.rows if row.label not in set(taxa)]
            return Alignment(keep, alignment.nchar, alignment.tail)

        self.addCleanup(STEPS.pop, "drop_taxa", None)
        STEPS["drop_taxa"] = drop_taxa

    def test_custom_step_runs_and_updates_ntax(self):
        alignment = read_alignment(self.sample_file())
        result = apply_steps(
            alignment,
            ["truncate", {"action": "drop_taxa", "taxa": ["t1"]}],
            length=5,
        )
        lines = format_alignment(result).split("\n")
        self.assertEqual(lines[0], "2  5")  # ntax recomputed
        self.assertEqual(lines[1], "OUT     ACGTA")
        self.assertEqual(lines[2], "t2      GGGGA")

    def test_unknown_step_is_reported(self):
        alignment = read_alignment(self.sample_file())
        with self.assertRaisesRegex(ValueError, "unknown step"):
            apply_steps(alignment, ["nope"], length=5)


class ConfigTests(TempDirTestCase):
    def test_resolves_relative_paths_and_sorts_lengths(self):
        path = self.tmp / "run.toml"
        path.write_text(
            'input_dir = "in"\n'
            'output_dir = "out"\n'
            "lengths = [200, 100, 100]\n"
        )
        config = load_config(path)
        self.assertEqual(config["input_dir"], (self.tmp / "in").resolve())
        self.assertEqual(config["output_dir"], (self.tmp / "out").resolve())
        self.assertEqual(config["lengths"], [100, 200])
        self.assertEqual(config["steps"], ["truncate"])  # default

    def test_reads_steps(self):
        path = self.tmp / "run.toml"
        path.write_text(
            'input_dir = "in"\n'
            'output_dir = "out"\n'
            "lengths = [100]\n"
            'steps = ["truncate", {action = "drop_taxa", taxa = ["t1"]}]\n'
        )
        self.assertEqual(
            load_config(path)["steps"],
            ["truncate", {"action": "drop_taxa", "taxa": ["t1"]}],
        )

    def test_rejects_a_config_without_lengths(self):
        path = self.tmp / "bad.toml"
        path.write_text('input_dir = "in"\noutput_dir = "out"\n')
        with self.assertRaises(ValueError):
            load_config(path)


class BuildDatasetsTests(TempDirTestCase):
    def test_writes_one_directory_per_length_and_mirrors_input_tree(self):
        seq_dir = self.tmp / "in" / "n25" / "03" / "seqs-500"
        seq_dir.mkdir(parents=True)
        (seq_dir / "0002.phy").write_text(SAMPLE)

        build_datasets(self.tmp / "in", self.tmp / "out", [4, 8])

        for length in (4, 8):
            target = self.tmp / "out" / str(length) / "n25" / "03" / "seqs-500" / "0002.phy"
            self.assertTrue(target.is_file())
            header, first_row = target.read_text().split("\n")[:2]
            self.assertEqual(header, f"3  {length}")
            self.assertEqual(first_row, "OUT     " + "ACGTACGTACGT"[:length])

    def test_skips_files_it_cannot_parse(self):
        seq_dir = self.tmp / "in"
        seq_dir.mkdir()
        (seq_dir / "good.phy").write_text(SAMPLE)
        (seq_dir / "bad.phy").write_text("junk\n")

        build_datasets(self.tmp / "in", self.tmp / "out", [4])

        self.assertTrue((self.tmp / "out" / "4" / "good.phy").is_file())
        self.assertFalse((self.tmp / "out" / "4" / "bad.phy").exists())


if __name__ == "__main__":
    unittest.main()
