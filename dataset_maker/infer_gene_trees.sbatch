#!/usr/bin/env bash
# Infer IQ-TREE gene trees for every replicate of a dataset directory.
#
# For each seqs-500 folder under the input directory, every .phy alignment is
# given to IQ-TREE (-m GTR+G, one alignment per process, $SLURM_CPUS_PER_TASK
# processes running concurrently) and the ML trees are collected into a
# iqtree_500.nwk file next to the seqs-500 folder. Intermediates go into an
# iqtree-work/ folder next to each seqs-500 folder. Finished genes are
# skipped on rerun, so after a timeout or failure the job can simply be
# resubmitted as-is.
#
# Submit from the repository root (paths are resolved relative to it):
#   mkdir -p logs
#   sbatch dataset_maker/infer_gene_trees.sbatch
# or with a different input directory and extra infer_gene_trees.py options:
#   sbatch dataset_maker/infer_gene_trees.sbatch data/camus-dataset/n50 --clean
#
# Sizing: 50 replicates x 1000 genes per length, three lengths in
# data/trimmed, i.e. 150k short IQ-TREE runs; at ~5 s per gene on 32 cores
# this fits comfortably in one day.
#SBATCH --time=7-00:00:00
#SBATCH --job-name=infer-gene-trees
#SBATCH --partition=cs
#SBATCH --account=cs-warnow
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32g
#SBATCH --output=logs/infer_gene_trees_%j_log.txt

set -u

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

# First argument is the input directory (default: the trimmed datasets);
# any further arguments are passed through to infer_gene_trees.py.
INPUT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project/data/trimmed
shift 2>/dev/null || true

cd "$PROJECT"
mkdir -p logs

echo "=== input: $INPUT, threads: $SLURM_CPUS_PER_TASK ==="
"$PYTHON" dataset_maker/infer_gene_trees.py \
    --input "$INPUT" \
    --threads "$SLURM_CPUS_PER_TASK" \
    "$@"

echo "=== done ==="
