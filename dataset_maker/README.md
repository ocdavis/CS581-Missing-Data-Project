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

Optionally, it also produces missing-taxa variants: for each percentage `P`
in `missing_percents`, a full-length copy of the dataset is written under
`output_dir/miss_<P>/` in which every file has `P`% of its ingroup taxa
removed (rows dropped, header taxon count updated). The subset is drawn
independently for each file, the outgroup (`OUT` by default) is never removed
and is not counted toward the percentage, and the draw is seeded per file from
`seed` and the file's relative path, so the same seed always gives the same
dataset regardless of worker count. With 50 ingroup taxa, 20% removes 10 taxa
(41 rows left) and 50% removes 25 (26 rows left).

The two kinds of variant are not crossed: `len_*` copies keep every taxon and
`miss_*` copies keep the full sequence length.

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

# also write the missing-taxa variants
python3 dataset_maker/trim_phy.py \
    --input data/camus-dataset/n50 \
    --output data/trimmed \
    --lengths 100 250 500 \
    --missing 20 50 --outgroup OUT --seed 581
```

To produce only the missing-taxa variants, use a config with `lengths = []`.

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
non-zero if any replicate differs. On the 50 `n50` networks distributed with
BROOQS-Data, the computed ToBs match their `true_tob.nwk` files in all 50
replicates.

`tree_of_blobs(newick) -> newick` and
`bipartitions(newick) -> (leaves, splits)` can be imported for reuse, e.g.
when computing FN/FP rates. No third-party packages are required.
