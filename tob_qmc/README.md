# tob_qmc

Estimates trees of blobs with TOB-QMC (the tree-of-blobs mode of
[TREE-QMC](https://github.com/molloy-lab/TREE-QMC)) in two stages, each a
SLURM script that also runs with plain `bash`:

1. `run_base_tree.sbatch` estimates a TREE-QMC species tree from each
   gene-tree file. This base tree is the shared starting point of TOB-QMC
   and ECToBlob.
2. `run_tob_qmc.sbatch` tests every branch of the base tree for
   network-like signal and contracts the branches that fail, giving the
   tree of blobs.

Both work on a directory laid out as `<dir>/.../<rep>/iqtree_500.nwk` (all
gene trees of one replicate, one Newick tree per line) and write their
results next to the gene trees. FN/FP against the true trees of blobs are
computed with `../evaluation/compare_tob.py`.

## Base trees (`run_base_tree.sbatch`)

Estimates a TREE-QMC species tree for every gene-tree file under a
directory. This is the base ("refinement") tree that both TOB-QMC and
ECToBlob start from, so that the two methods contract the same tree. It is
`ectoblob/`'s `start_tree_glob = "tree-qmc.nwk"`.

For each `<dir>/iqtree_500.nwk` it runs

```bash
tree-qmc --root OUT -i <dir>/iqtree_500.nwk -o <dir>/tree-qmc.nwk --override
```

(TREE-QMC defaults: unweighted quartets, `n2` normalization) and writes
`tree-qmc.nwk` next to the gene trees, with TREE-QMC's messages in
`tree-qmc.log`. The trees are binary, contain every taxon of the replicate
and are rooted at `OUT`. Existing trees are skipped, failures are listed, and
the exit status is non-zero unless every gene-tree file has a tree.

It then adds branch lengths to each tree without changing its topology
(`tree-qmc --supportonly tree-qmc.nwk -i iqtree_500.nwk`):

| File | Contents |
|------|----------|
| `tree-qmc.nwk` | topology only |
| `tree-qmc.support.nwk` | TREE-QMC's annotated tree: quartet support labels (`'[f1=...;q1=...;EN=...]'`) and branch lengths |
| `tree-qmc.brlen.nwk` | the same tree with the labels removed, i.e. plain Newick with branch lengths |

The lengths are in coalescent units, estimated from the quartet frequencies
around each branch. Only internal branches have one; terminal branches have
no length. A length of `0` means the gene trees do not favour that branch
over the two alternative resolutions. Running the script again on a
directory that already has `tree-qmc.nwk` only adds the missing
branch-length files.

### Build

Base trees do not need tree-of-blobs support, so TREE-QMC can be built
without R. The current TREE-QMC `main` does not compile this way; commit
`7f840e7` (version 4.1.5) does:

```bash
module load gcc/13.3.0
git clone https://github.com/molloy-lab/TREE-QMC    # next to this repository
cd TREE-QMC && git checkout 7f840e7
cd external/MQLib && make
cd ../../ && mkdir -p build && cd build
g++ -std=c++11 -O2 -I ../external/MQLib/include -I ../external/toms743 \
    -o tree-qmc ../src/*.cpp ../external/toms743/toms743.cpp \
    ../external/MQLib/bin/MQLib.a -lm -DVERSION=\"$(cat ../version.txt)\"
./tree-qmc -h | head -5
```

The script looks for `../TREE-QMC/build/tree-qmc` relative to this
repository; set `TREE_QMC_BUILD` to a different `build` directory if needed.

### Run

With the gene trees laid out as `data/gene_trees/<dataset>/<rep>/iqtree_500.nwk`,
from the repository root:

```bash
# one replicate, without SLURM (about 7 s)
bash tob_qmc/run_base_tree.sbatch data/gene_trees/miss_50_len_500/00

# everything (250 trees, a few minutes on 16 cores)
mkdir -p logs
sbatch tob_qmc/run_base_tree.sbatch data/gene_trees
```

### Other datasets

The script works on any directory: it finds every gene-tree file below it,
at any depth, and writes the species tree next to each one. Each gene-tree
file must hold all gene trees of one dataset, one Newick tree per line. The
file names and the outgroup are environment variables:

| Variable | Default | Meaning |
|----------|---------|---------|
| `GENES_NAME` | `iqtree_500.nwk` | name of the gene-tree files to look for (a `find -name` pattern, so `"*.tre"` also works) |
| `OUT_NAME` | `tree-qmc.nwk` | name of the species tree written next to each one |
| `OUTGROUP` | `OUT` | taxon to root at; set `OUTGROUP=` (empty) to leave the tree unrooted |

```bash
GENES_NAME=genes.tre OUT_NAME=species.tre OUTGROUP=Homo_sapiens \
    sbatch tob_qmc/run_base_tree.sbatch path/to/dataset
```

With a pattern that matches several files in one folder, they would all
write to the same `OUT_NAME`, so keep one gene-tree file per folder.

To collect only the base trees, keeping the layout:

```bash
cd data && find gene_trees -name "tree-qmc*.nwk" | sort | tar -czf ../tree_qmc_trees.tar.gz -T -
```

## Trees of blobs (`run_tob_qmc.sbatch`)

For each folder holding `iqtree_500.nwk` and its base tree `tree-qmc.nwk`:

```bash
# 1. search for the minimum p-value around every branch of the base tree
tree-qmc --blobsearchonly tree-qmc.nwk --store_pvalue --iter_limit_blob 0 \
    -i iqtree_500.nwk -o tob-qmc.pvalues.nwk

# 2. contract branches at the chosen significance levels
tree-qmc --alpha 1e-7 --beta 0.95 --load_pvalue \
    -i tob-qmc.pvalues.nwk -o tob-qmc.nwk
```

| File | Contents |
|------|----------|
| `tob-qmc.nwk` | estimated tree of blobs (blobs are polytomies) |
| `tob-qmc.pvalues.nwk` | base tree with, on every branch, the tree-test p-value (`qtt_p`), star-test p-value (`qst_p`), quartet counts and the four taxa giving the minimum |
| `tob-qmc.search.log`, `tob-qmc.contract.log` | tree-qmc's messages for the two steps |

`--iter_limit_blob 0` is the exhaustive search (every 4-taxon set across a
branch); alpha and beta are TOB-QMC's defaults. This is the three-step form
of `tree-qmc --blob` from the TREE-QMC tutorial, with the base tree taken
from `run_base_tree.sbatch` instead of being re-estimated.

Finished steps are skipped, so an interrupted job can simply be
resubmitted, and failures are listed without stopping the rest. The search
is the slow step: about 9 minutes per replicate for 51 taxa and 1000 gene
trees, so roughly 1 to 1.5 hours for 250 replicates on 32 cores.

### Build

This stage needs TREE-QMC built **with** tree-of-blobs support, which
requires cmake 3.18+, a C++20 compiler and R with Rcpp, RInside and
MSCquartets. On the campus cluster the system gcc, cmake and R are enough;
do this with no modules loaded (`module purge`), so that R, the R packages
and tree-qmc all use the same compiler.

```bash
# R packages, into R's personal library (slow: 20-40 minutes)
Rscript -e "lib <- Sys.getenv('R_LIBS_USER'); dir.create(lib, recursive=TRUE, showWarnings=FALSE); install.packages(c('Rcpp','RInside','MSCquartets'), lib=lib, repos='https://cloud.r-project.org/')"
Rscript --vanilla -e "for (p in c('Rcpp','RInside','MSCquartets')) cat(p, requireNamespace(p, quietly=TRUE), '\n')"

# TREE-QMC, in its own folder next to this repository
git clone https://github.com/molloy-lab/TREE-QMC TREE-QMC-tob
cd TREE-QMC-tob && git checkout 53b15e5
mkdir -p build && cd build
cmake .. && cmake --build . -j 4
./tree-qmc -h | head -5
```

Commit `53b15e5` (22 Sep 2026, version 4.1.5) is the one used for this
project. The script looks for `../TREE-QMC-tob/build/tree-qmc` relative to
this repository; set `TREE_QMC_BUILD` to another `build` directory if
needed. The R library must not be set through `~/.Renviron`, which the
build ignores.

### Run

From the repository root, with no modules loaded:

```bash
# one replicate, without SLURM (about 9 minutes)
bash tob_qmc/run_tob_qmc.sbatch data/gene_trees/miss_0_len_500/00

# everything
mkdir -p logs
sbatch --partition=cs,secondary --time=04:00:00 tob_qmc/run_tob_qmc.sbatch data/gene_trees
```

As with the base trees, the `#SBATCH` lines use the `cs-melkebir` account;
pass `--account=NAME` before the script name to change it.

### Settings

Environment variables, set before `sbatch` or `bash`:

| Variable | Default | Meaning |
|----------|---------|---------|
| `ALPHA` | `1e-7` | significance level of the tree test |
| `BETA` | `0.95` | significance level of the star test |
| `SEARCH` | `--iter_limit_blob 0` | search options; `SEARCH=` (empty) for TOB-QMC's default heuristic, `SEARCH=--3f1a` for 3-fix-1-alter |
| `GENES_NAME` | `iqtree_500.nwk` | gene-tree file name |
| `BASE_NAME` | `tree-qmc.nwk` | base tree file name |
| `PVAL_NAME` | `tob-qmc.pvalues.nwk` | p-value tree file name |
| `TOB_NAME` | `tob-qmc.nwk` | tree of blobs file name |

Because the p-values are kept, other alpha/beta values only need the
contraction step, which takes seconds. Give the result a different name so
the default one is kept:

```bash
ALPHA=1e-6 TOB_NAME=tob-qmc.a1e-6.nwk bash tob_qmc/run_tob_qmc.sbatch data/gene_trees
```

To collect the results, keeping the layout:

```bash
cd data && find gene_trees -name "tob-qmc*.nwk" | sort | tar -czf ../tob_qmc_trees.tar.gz -T -
```
