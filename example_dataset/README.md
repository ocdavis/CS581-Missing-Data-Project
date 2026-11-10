# example_dataset

A small slice of the study dataset, small enough to keep in the repository.
It has exactly the layout of the full dataset and was produced with the same
two commands, just on fewer files: replicates `00` and `01` of the 50-taxon
(+ outgroup) condition of the
[CAMUS dataset](https://doi.org/10.13012/B2IDB-6892704_V1), with only the
first 10 of the 1000 genes of each replicate.

## Layout

```text
example_dataset/
├── miss_0_len_100/<rep>/seqs-500/*.phy    all 51 taxa
├── miss_0_len_200/ ...
├── miss_0_len_500/ ...                    the original alignments
├── miss_20_len_100/ ...                   20% of ingroup taxa removed per gene
├── miss_20_len_200/ ...                   (41 taxa left)
├── miss_20_len_500/ ...
├── miss_50_len_100/ ...                   50% of ingroup taxa removed per gene
├── miss_50_len_200/ ...                   (26 taxa left)
├── miss_50_len_500/ ...
└── true_tob/<rep>/true_tob.nwk            true tree of blobs
```

`miss_<P>_len_<L>` holds the alignments with `P`% of the ingroup taxa removed
and the sequences truncated to `L` sites. `miss_0_len_500` is the original
data (same taxa and sequences as CAMUS, rewritten one taxon per line).

- Taxa are removed independently for each gene, so every taxon still occurs
  in each replicate. `true_tob/<rep>/true_tob.nwk` is therefore the reference
  for **every** condition of that replicate.
- For a given gene and percentage, the same taxa are removed at every length.
- The outgroup `OUT` is never removed.

## Generating the full dataset

Download `camus-dataset.tar.xz` from the CAMUS dataset page and extract the
`n50` condition so that `data/camus-dataset/n50/00/seqs-500/0001.phy` exists:

```bash
mkdir -p data
tar -xJf camus-dataset.tar.xz -C data camus-dataset/n50
```

Then, from the repository root (Python 3.11+, no packages needed):

```bash
python dataset_maker/trim_phy.py --config dataset_maker/example_config.toml
python dataset_maker/true_tob.py --input data/camus-dataset/n50 --output data/trimmed/true_tob
```

This writes the layout above under `data/trimmed/` for all 50 replicates and
1000 genes (about 6 GB; `data/` is gitignored).

The taxon removal is seeded (`seed = 581` in the config) and depends only on
the seed and each file's path below the input directory, so the full run
produces, for these 20 genes, exactly the files stored here. Everyone who
runs it gets the same dataset.

## How this example was made

The same two commands, with a copy of the first 10 alignments (and
`true_net.nwk`) of replicates `00` and `01` as the input directory:

```bash
python dataset_maker/trim_phy.py --input <subset> --output example_dataset \
    --lengths 100 200 500 --missing 0 20 50 --outgroup OUT --seed 581
python dataset_maker/true_tob.py --input <subset> --output example_dataset/true_tob
```
