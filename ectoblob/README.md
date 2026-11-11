# ectoblob

Infers a tree of blobs (ToB) with **ECToBlob** ([MSCquartets R package](
https://cran.r-project.org/web/packages/MSCquartets/vignettes/ECToBlob.html))
for every gene-tree file in a directory — the ECToBlob counterpart of
`../tob_qmc`.

ECToBlob needs two inputs per replicate:

1. **Gene trees** — a multi-tree Newick file, one tree per line (e.g. the
   `iqtree_500.nwk` files produced by `../dataset_maker/infer_gene_trees.py`).
   Missing taxa are allowed, but every set of 4 taxa must appear in at least
   one gene tree.
2. **A starting tree** assumed to be a resolution of the true tree of blobs —
   e.g. a TREE-QMC species tree (run *without* `--blob`) or an ASTRAL tree.
   The starting tree must not contain polytomies, so a TOB-QMC output is not
   suitable. For exploratory purposes only, ECToBlob's own
   `quartetDistTree()` can build one from the gene trees.

ECToBlob tests each quartet for fit to a star tree and to the resolution
displayed on the starting tree, contracts unsupported edges, and returns a
sequence of trees; the script selects the `indexEarly` (or `indexLate`) tree
from that sequence and roots it at the outgroup.

## `run_ectoblob.py`

For each file matching `gene_tree_glob` under `input_dir`, the starting tree
is looked up with `start_tree_glob` in the same directory (exactly one match
per replicate is required) and `ectoblob.R` is invoked via `Rscript`. The
selected ToB is written to `<output_dir>/<relative path>/<stem>.ectoblob.tre`
(with R's messages in `<stem>.ectoblob.log` next to it); with
`output_dir = input_dir` the tree lands next to the gene trees. Files whose
output tree already exists are skipped, and the exit status is non-zero if
any replicate failed.

### Usage

```bash
python3 ectoblob/run_ectoblob.py --config ectoblob/example_config.toml
```

### Configuration

See `example_config.toml` for all keys. The analysis parameters follow the
ECToBlob vignette: `alpha` (which tree in the sequence is selected), `beta`
(which edges are collapsed as unresolved up front), `qtype` (`bi`, `quad`,
or `mul`), `test_correction` (`Bon`, `Cauchy`, `BBC`, or `CBC`), and `index`
(`early` or `late`). A reported analysis should state all of these, and
exploring multiple settings is recommended — set `save_table = true` to keep
the quartet p-value table (`<stem>.pTable.rds`) so later runs with different
parameters skip the slow quartet tallying.

### Prerequisites

R with the `MSCquartets` package (and its dependencies `ape` and
`phangorn`); `rscript` names the R executable (default `Rscript` on PATH).
