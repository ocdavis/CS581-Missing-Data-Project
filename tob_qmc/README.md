# tob_qmc

Runs TOB-QMC (the `--blob` mode of
[TREE-QMC](https://github.com/molloy-lab/TREE-QMC)) on every gene-tree file
in a directory, with TOB-QMC's default alpha (1e-7) and beta (0.95) and the
exhaustive p-value search.

## 1. Build tree-qmc

TOB-QMC needs a `tree-qmc` binary built **with tree-of-blobs support**. This
is a Linux build (use the cluster, or WSL on Windows). It needs git, cmake
3.18+, a C++ compiler and R.

Install the R packages (takes a while):

```bash
Rscript -e "install.packages(c('Rcpp','RInside','MSCquartets'), repos='https://cloud.r-project.org/')"
```

Build:

```bash
git clone https://github.com/molloy-lab/TREE-QMC.git
cd TREE-QMC
mkdir -p build && cd build
cmake ..
cmake --build . -j
```

The binary is `TREE-QMC/build/tree-qmc`. Either add that `build` directory to
your `PATH`, or put the full path in `tree_qmc` in the config. More detail is
in the
[TREE-QMC tree-of-blobs tutorial](https://github.com/molloy-lab/TREE-QMC/blob/main/tutorial/tree-of-blobs/README.md).

## 2. Edit the config

Copy `example_config.toml` and set:

- `input_dir` and `gene_tree_glob`: where the gene-tree files are. Each
  matching file is one input (all gene trees of one replicate, one Newick
  tree per line). The example points at the true gene trees in the CAMUS
  dataset (`data/camus-dataset/n50/*/g_true.nwk`).
- `output_dir`: where the estimated ToBs go.

## 3. Run

From the repository root (needs Python 3.11+, no packages):

```bash
python3 tob_qmc/run_tob_qmc.py --config tob_qmc/example_config.toml
```

For each gene-tree file `<input_dir>/<rep>/<name>.nwk` this runs

```bash
tree-qmc --blob --alpha 1e-07 --beta 0.95 --root OUT \
    -i <input_dir>/<rep>/<name>.nwk \
    -o <output_dir>/<rep>/<name>.tob.tre --override --iter_limit_blob 0
```

and writes:

| File | Contents |
|------|----------|
| `<output_dir>/<rep>/<name>.tob.tre` | estimated tree of blobs |
| `<output_dir>/<rep>/<name>.tob.log` | tree-qmc's messages for that run |

Files whose `.tob.tre` already exists are skipped, so the same command can be
rerun after an interruption. Delete an output to redo it. If `tree-qmc` fails
on a file, the script continues with the rest and reports the failures at the
end.

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
