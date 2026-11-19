# evaluation

## `compare_tob.py`

Counts false negative (FN) and false positive (FP) branches of estimated
trees against reference trees. It works for any pair of Newick trees: a tree
of blobs from TOB-QMC or ECToBlob against the true tree of blobs, one method
against another, or trees from a different dataset.

Trees are compared as sets of non-trivial unrooted bipartitions (internal
branches):

| | |
|---|---|
| FN | branches of the reference that are missing from the estimate |
| FP | branches of the estimate that are not in the reference |
| `fn_rate` | FN / internal branches of the reference |
| `fp_rate` | FP / internal branches of the estimate (0 if it has none) |

Rooting, branch lengths, support values, internal node labels, `[comments]`
and quoted annotations (such as TREE-QMC's `'[f1=...]'` and TOB-QMC's
`'[blob_id=...]'`) are ignored, and polytomies are allowed. Only the first
tree of each file is read. Needs Python 3.9+ and no packages.

### Two trees

```bash
python3 evaluation/compare_tob.py --estimated estimate.tre --reference reference.tre
```

### A whole directory

Every file under `--estimated` whose name matches one of the `--glob`
patterns is compared with its reference. The reference is looked up under
`--reference`, first in the same relative folder as the estimate and then in
a folder named like the estimate's own folder (the replicate). Its file name
is `--reference-name` (default `true_tob.nwk`).

For this project, with the true trees of blobs made by
`dataset_maker/true_tob.py --input data/camus-dataset/n50`
(`data/camus-dataset/n50/<rep>/true_tob.nwk`):

```bash
# TOB-QMC (data/gene_trees/<condition>/<rep>/tob-qmc.nwk)
python3 evaluation/compare_tob.py --estimated data/gene_trees --glob tob-qmc.nwk \
    --reference data/camus-dataset/n50 --output results/tob_qmc_errors.csv

# ECToBlob (results/ectoblob/<condition>/<rep>/iqtree_500.ectoblob.tre)
python3 evaluation/compare_tob.py --estimated results/ectoblob --glob "*.ectoblob.tre" \
    --reference data/camus-dataset/n50 --output results/ectoblob_errors.csv

# several kinds of tree at once, e.g. the base trees as a baseline
python3 evaluation/compare_tob.py --estimated data/gene_trees \
    --glob tob-qmc.nwk tree-qmc.nwk "*.ectoblob.tre" \
    --reference data/camus-dataset/n50 --output results/errors.csv
```

Quote patterns containing `*` so the shell does not expand them.

### One method against another

Point `--reference` at the other method's trees. With both in the same
folders:

```bash
python3 evaluation/compare_tob.py --estimated data/gene_trees --glob "*.ectoblob.tre" \
    --reference data/gene_trees --reference-name tob-qmc.nwk
```

Here "FN" are branches TOB-QMC kept and ECToBlob contracted, and "FP" the
reverse. This is also the way to compare methods on a dataset without a known
true tree (e.g. an empirical one), where FN/FP are relative to whichever tree
is given as the reference rather than errors.

### Different taxon sets

By default the two trees must have exactly the same taxa, and a mismatch is
reported as an error for that file. With `--common-taxa`, both trees are
first restricted to the taxa they share (at least 4), and `n_taxa` in the
output is the number of shared taxa.

### Output

One CSV row per estimated tree, to `--output` or to the screen:

| Column | Meaning |
|---|---|
| `method` | file name of the estimated tree |
| `condition` | folders above the replicate folder (e.g. `miss_50_len_500`) |
| `replicate` | folder containing the estimated tree |
| `n_taxa` | taxa compared |
| `ref_edges`, `est_edges` | internal branches of the reference and the estimate |
| `fn`, `fp`, `fn_rate`, `fp_rate` | as defined above |
| `file`, `reference` | the two files compared |

A table of mean FN, FP and rates per method and condition is also printed.
Files that cannot be compared (missing reference, unreadable tree, different
taxa) are listed, the rest are still written, and the exit status is
non-zero.

## `polytomy_stats.py`

Number of polytomies and largest polytomy degree of each estimated tree and
of its reference, found with the same options as `compare_tob.py`. A polytomy
is an internal node of degree 4 or more in the unrooted tree (one blob); its
degree is the number of subtrees attached to it.

```bash
python3 evaluation/polytomy_stats.py --estimated data/gene_trees \
    --glob tob-qmc.nwk "*.ectoblob.tre" \
    --reference data/camus-dataset/n50 --output results/polytomy_stats.csv
```

## `plot_errors.py`

Box plots of FN and FP per condition, one box per method, from one or
more CSV files written by `compare_tob.py` (needs matplotlib):

```bash
python3 evaluation/plot_errors.py results/base_tree_errors.csv \
    results/tob_qmc_errors.csv results/ectoblob_errors.csv \
    --output results/fn_fp.png --title "CAMUS n50"
```

Rates are plotted by default and `--counts` plots branches per tree.
`--conditions miss_0_len_500 miss_0_len_200 miss_0_len_100` limits the
figure to the named conditions, `--xlabel` sets the text under the x axis,
and the extension of `--output` picks the format (`.png`, `.pdf`, `.svg`).
Both panels share one scale. A method with no rows for a condition has no
box there.

## `plot_overview.py`

All experiments on one page: one row per experiment (sequence length, taxa
missing from every gene, 50% and 80% of taxa missing from a share of the
genes), FN rate on the left and FP rate on the right, on one scale.

```bash
python3 evaluation/plot_overview.py results/base_tree_errors.csv \
    results/tob_qmc_errors.csv results/ectoblob_errors.csv \
    --output results/overview.png
```

The scale stops at 70% by default; values above it (trees with no internal
branches have a false negative rate of 100%) are drawn as triangles at the
top edge. `--top 100` shows the full scale. The rows and their conditions
are listed at the top of the script.

## Result files (`../results`)

| File | Contents |
|---|---|
| `tob_qmc_errors.csv` | FN/FP of TOB-QMC, 16 conditions x 50 replicates |
| `base_tree_errors.csv` | FN/FP of the TREE-QMC base trees, same conditions |
| `ectoblob_errors.csv` | FN/FP of ECToBlob, the 13 conditions it ran on |
| `tob_qmc_polytomies.csv` | polytomy counts and degrees of TOB-QMC and the true trees |
| `overview.png`, `fn_fp_*.png` | figures made by the two plotting scripts |

The base-tree rows of the nine conditions added last (25 and 50 sites, 95%
and 99% of genes, 80% of taxa missing) were computed from
`tob-qmc.pvalues.nwk`, which is the base tree with p-value labels.

Contracting can only remove branches, so for a tree of blobs built from a
base tree the base tree's FN is the lowest FN the method can reach; running
the script on the base trees (`tree-qmc.nwk`) gives that floor.
