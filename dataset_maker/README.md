# dataset_maker

Tools for generating trimmed variants of the CAMUS PHYLIP datasets.

## `trim_phy.py`

Reads every `.phy` file under an input directory, truncates each sequence to a
target length, and writes the result to an output directory while preserving
the relative file structure. One full copy of the dataset is produced per
length, under `output_dir/len_<L>/`.

Handles both sequential and interleaved PHYLIP layouts. Sequences shorter
than the target length are left unchanged (no padding).

Files are processed in parallel with a process pool (`--workers`, default:
all cores); each file is parsed once and written once per length. A progress
bar is shown when tqdm is installed, otherwise a plain counter is printed.

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
```

### Extending

The pipeline is deliberately modular — each stage is a plain function you can
import and reuse or replace:

| Stage          | Function                          |
|----------------|-----------------------------------|
| parse          | `parse_phylip(path) -> PhylipFile` |
| transform      | `trim_sequences(data, length) -> PhylipFile` |
| write          | `write_phylip(data, path)`         |
| per-file driver| `process_file(src, input_dir, output_dir, lengths)` |
| whole-run driver| `run(input_dir, output_dir, lengths, workers)` |

`PhylipFile` is a small dataclass holding `n_taxa`, `n_sites`, and
`rows: list[tuple[name, sequence]]` — add fields or new transform functions
alongside `trim_sequences` to add your own improvements.
