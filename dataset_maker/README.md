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
(percentage `P`, length `L`) pair under `output_dir/miss_<P>/len_<L>/`, in
which every file has `P`% of its ingroup taxa removed (rows dropped, header
taxon count updated) and its sequences truncated to `L`. For example,
`missing_percents = [0, 20, 50]` with `lengths = [100, 200, 500]` gives nine
datasets, from `miss_0/len_100/` to `miss_50/len_500/`. Include `0` to get
copies with all taxa, and the original length (500) to get untrimmed copies;
`miss_0/len_500/` is then a copy of the input alignments, so every condition
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

# cross the lengths with missing-taxa percentages (miss_<P>/len_<L>/)
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
