#!/usr/bin/env bash
# Generate the two partial-missingness datasets at full sequence length:
#
#   data/trimmed/miss_50_files_50_len_500   50% of the genes lose 50% of
#                                           their ingroup taxa, the rest
#                                           keep all taxa
#   data/trimmed/miss_50_files_80_len_500   80% of the genes lose 50% of
#                                           their ingroup taxa, the rest
#                                           keep all taxa
#
# Every sequence stays 500 sites long. All logic and every parameter lives in
# dataset_maker/partial_miss_config.toml; this wrapper only sets the SLURM
# resources and points trim_phy.py at that config. Existing datasets under
# data/trimmed are untouched: only the two new condition folders are written.
#
# This step is pure file I/O (50k alignments read, 100k written, ~1.6 GB), so
# it needs far fewer resources than the gene-tree jobs. To infer gene trees
# afterwards, run e.g.
#   sbatch dataset_maker/infer_gene_trees.sbatch \
#       data/trimmed/miss_50_files_50_len_500
#
# Submit from the repository root:
#   sbatch dataset_maker/run_partial_miss.sbatch
#SBATCH --time=04:00:00
#SBATCH --job-name=trim-partial-miss
#SBATCH --partition=secondary
#SBATCH --account=campusclusterusers
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g
#SBATCH --output=logs/trim_partial_miss_%j_log.txt

set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

cd "$PROJECT"
mkdir -p logs

"$PYTHON" dataset_maker/trim_phy.py \
    --config dataset_maker/partial_miss_config.toml \
    --workers "$SLURM_CPUS_PER_TASK"

# Informational only: trim_phy already failed above if anything went wrong.
echo "=== datasets written under data/trimmed ==="
ls -d data/trimmed/miss_50_files_*_len_500 2>/dev/null || true
