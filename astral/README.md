# astral

Estimates an ASTRAL-III species tree for every gene-tree file under a
directory, with ASTRAL's default options:

```bash
java -jar astral.5.7.8.jar -i <dir>/iqtree_500.nwk -o <dir>/astral.nwk
```

The jar is `repos/ASTRAL/Astral/astral.5.7.8.jar`; ASTRAL's messages go to
`<dir>/astral.log`. Trees are written next to the gene trees, so the layout
matches the other methods:

```
data/gene_trees/<dataset>/<rep>/iqtree_500.nwk
data/gene_trees/<dataset>/<rep>/astral.nwk
```

To gather them in a separate directory afterwards, keeping the subdirectory
structure, run `collect_astral_trees.sh` (see below).

Existing trees are skipped and failures are listed, so an interrupted run can
simply be restarted and the exit status is non-zero unless every gene-tree
file has a tree.

## Run

One SLURM job per dataset directory, submitted from the repository root:

```bash
bash astral/run_all_astral.sh                          # all 16 datasets
bash astral/run_all_astral.sh miss_0_len_500 miss_50_len_500
```

`run_all_astral.sh` skips directories that already have every `astral.nwk`,
so relaunching it after an interruption only resubmits the missing work.
Each job runs the replicates of its dataset one after another (ASTRAL is
single-threaded, one core is enough).

Or run a single directory, with or without SLURM:

```bash
sbatch astral/astral_one.sbatch data/gene_trees/miss_0_len_500
bash astral/astral_one.sbatch data/gene_trees/miss_0_len_500
python3 astral/run_astral.py data/gene_trees/miss_0_len_500/00
```

## Collect the trees

`collect_astral_trees.sh` copies the trees written next to the gene trees
into `data/astral_trees`, keeping the subdirectory structure:

```bash
bash astral/collect_astral_trees.sh
# data/gene_trees/<dataset>/<rep>/astral.nwk
#   -> data/astral_trees/<dataset>/<rep>/astral.nwk
```

It can run at any point, so it also works as a snapshot of the trees that
exist so far; pass two arguments to use other directories.

## Other datasets

`run_astral.py` searches any directory recursively, at any depth. The file
names are options (or environment variables for the driver script):

| Option | Driver variable | Default | Meaning |
|--------|-----------------|---------|---------|
| `--pattern` | `GENES_NAME` | `iqtree_500.nwk` | gene-tree file name to look for |
| `--out-name` | `OUT_NAME` | `astral.nwk` | species-tree file name |
| `--jar` | | `repos/ASTRAL/Astral/astral.5.7.8.jar` | ASTRAL jar |

Use `--output DIR` to write the species trees somewhere other than next to
the gene trees, keeping the path relative to the input directory:

```bash
python3 astral/run_astral.py data/gene_trees --output data/astral_trees
# data/gene_trees/miss_0_len_500/00/iqtree_500.nwk
#   -> data/astral_trees/miss_0_len_500/00/astral.nwk
```
