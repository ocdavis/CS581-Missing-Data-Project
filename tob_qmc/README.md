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
