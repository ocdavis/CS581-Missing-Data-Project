# dataset_maker

Tools for generating trimmed and missing-taxa variants of the CAMUS PHYLIP
datasets.

## `trim_phy.py`

Reads every `.phy` file under an input directory, truncates each sequence to a
target length, and writes the result to an output directory while preserving
the relative file structure. One full copy of the dataset is produced per
length, under `output_dir/len_<L>/`.

Handles both sequential and interleaved PHYLIP layouts. Sequences shorter
than the target length are left unchanged (no padding).

Optionally, lengths are crossed with missing-taxa percentages: when
`missing_percents` is given, one copy of the dataset is written for every
(percentage `P`, length `L`) pair under `output_dir/miss_<P>_len_<L>/`, in
which every file has `P`% of its ingroup taxa removed (rows dropped, header
taxon count updated) and its sequences truncated to `L`. For example,
`missing_percents = [0, 20, 50]` with `lengths = [100, 200, 500]` gives nine
datasets, from `miss_0_len_100/` to `miss_50_len_500/`. Include `0` to get
copies with all taxa, and the original length (500) to get untrimmed copies;
`miss_0_len_500/` is then a copy of the input alignments, so every condition
lives in the same directory layout.

The subset of taxa is drawn independently for each file, and the same subset
is used at every length. The outgroup (`OUT` by default) is never removed and
is not counted toward the percentage. The draw is seeded per file from `seed`
and the file's relative path, so the same seed always gives the same dataset
regardless of worker count. With 50 ingroup taxa, 20% removes 10 taxa (41
rows left) and 50% removes 25 (26 rows left).

Without `missing_percents`, the output layout is `output_dir/len_<L>/` as
before.

Files are processed in parallel with a process pool (`--workers`, default:
all cores); each file is parsed once and written once per length. Files are
streamed from the directory walk with only a bounded window of tasks in
flight, so memory stays flat regardless of dataset size. A progress bar is
shown when tqdm is installed, otherwise a plain counter is printed.

### Usage

```bash
# via config file
python3 dataset_maker/trim_phy.py --config dataset_maker/example_config.toml

# or directly
python3 dataset_maker/trim_phy.py \
    --input data/camus-dataset/n50 \
    --output data/trimmed \
    --lengths 100 250 500 \
    --workers 8

# cross the lengths with missing-taxa percentages (miss_<P>_len_<L>/)
python3 dataset_maker/trim_phy.py \
    --input data/camus-dataset/n50 \
    --output data/trimmed \
    --lengths 100 250 500 \
    --missing 0 20 50 --outgroup OUT --seed 581
```
### Extending

The pipeline is deliberately modular — each stage is a plain function you can
import and reuse or replace:

| Stage          | Function                          |
|----------------|-----------------------------------|
| parse          | `parse_phylip(path) -> PhylipFile` |
| transform      | `trim_sequences(data, length) -> PhylipFile` |
| transform      | `remove_random_taxa(data, percent, rng, outgroup) -> PhylipFile` |
| write          | `write_phylip(data, path)`         |
| per-file driver| `process_file(src, input_dir, output_dir, lengths, missing, outgroup, seed)` |
| whole-run driver| `run(input_dir, output_dir, lengths, workers, missing, outgroup, seed)` |

`PhylipFile` is a small dataclass holding `n_taxa`, `n_sites`, and
`rows: list[tuple[name, sequence]]` — add fields or new transform functions
alongside `trim_sequences` to add your own improvements.

## `true_tob.py`

Computes the true tree of blobs (ToB) of every `true_net.nwk` under an input
directory and writes it as `true_tob.nwk` next to the network (or under
`--output`, mirroring the directory layout).

The ToB is obtained from the undirected graph of the network by finding its
cut edges (bridges), contracting every other edge so that each blob becomes a
single node, and suppressing the degree-2 nodes that remain. Blobs therefore
show up as polytomies. Branch lengths of the surviving edges are kept. The
tree is written rooted where the network was rooted, but should be treated as
unrooted when comparing.

Because taxa are removed per gene in the `miss_*` datasets, every taxon still
occurs in each replicate, so the same true ToB is the reference for all
missing-data and sequence-length conditions.

### Usage

```bash
python3 dataset_maker/true_tob.py --input data/camus-dataset/n50

# also compare against reference ToBs at the same relative paths, e.g. the
# true_tob.nwk files from https://github.com/shayesteh99/BROOQS-Data
python3 dataset_maker/true_tob.py --input data/camus-dataset/n50 \
    --reference path/to/BROOQS-Data/fixed-ILS-dataset/n50
```

The comparison is by sets of unrooted bipartitions; the exit status is
non-zero if any replicate differs.

`tree_of_blobs(newick) -> newick` and
`bipartitions(newick) -> (leaves, splits)` can be imported for reuse, e.g.
when computing FN/FP rates.

## `infer_gene_trees.py`

Infers gene trees with IQ-TREE for every replicate of a dataset directory.
Requires the `iqtree` binary (installed via `bioconda::iqtree` in
`environment.yml`).

Every `seqs-500` folder found under the input directory (e.g.
`data/trimmed/len_100/00/seqs-500`) is treated as one replicate: each `.phy`
alignment inside is given to IQ-TREE and the ML trees are collected into a
single multi-tree Newick file, one tree per line, written next to the
`seqs-500` folder (or under `--output`, mirroring the directory layout). The
output is named `g_500.nwk` by default, matching the gene-tree files
distributed with the CAMUS datasets, regardless of the actual alignment
length.

IQ-TREE runs one alignment at a time (`-T 1`) with `--threads` processes
running concurrently; for the many small alignments in these datasets this
is faster than multi-threaded runs. Model selection uses `-m MFP` unless
`--model` gives a fixed model, and `--seed` makes runs reproducible.

A gene is only (re-)run when its `<stem>.treefile` does not exist yet, so an
interrupted run can simply be restarted; pass `--overwrite` to redo
everything. Intermediate IQ-TREE files (`.log`, `.iqtree`, `.ckp.gz`, ...)
go into an `iqtree-work/` folder next to each `seqs-500` folder (or under
`--workdir DIR`, mirroring the input layout); `--clean` deletes each
replicate's work folder after it finishes successfully. A replicate only
gets its output file when every alignment in it produced a tree; otherwise
it is reported as failed and the exit status is non-zero.

### Usage

```bash
# via config file
python3 dataset_maker/infer_gene_trees.py \
    --config dataset_maker/infer_gene_trees_example.toml

# or directly
python3 dataset_maker/infer_gene_trees.py --input data/trimmed/len_100

# 16 concurrent IQ-TREE processes, fixed model, custom output location
python3 dataset_maker/infer_gene_trees.py --input data/camus-dataset/n50 \
    --threads 16 --model GTR+G --output data/gene_trees
```

### `infer_gene_trees.sbatch`

Submits the inference as a SLURM job (32 cores, one day limit), running the
script over the whole `data/trimmed` directory by default:

```bash
# from the repository root
mkdir -p logs
sbatch dataset_maker/infer_gene_trees.sbatch

# different input directory, plus extra infer_gene_trees.py options
sbatch dataset_maker/infer_gene_trees.sbatch data/camus-dataset/n50 --clean
```

The first argument is the input directory (default `data/trimmed`); further
arguments are passed through to `infer_gene_trees.py`. Because finished
genes are skipped, a timed-out or failed job can simply be resubmitted.
