#!/usr/bin/env bash
# Run ECToBlob on the collected gene trees as a SLURM job array with one
# task per dataset: task N processes the Nth dataset directory under
# data/gene_trees (miss_0_len_100, miss_0_len_200, ...). For each
# iqtree_500.nwk in its dataset, the task infers a tree of blobs with the
# TREE-QMC species tree next to it as the starting tree, writing
# <stem>.ectoblob.tre under results/ectoblob (mirroring the layout).
# Finished replicates are skipped, so the job can simply be resubmitted.
#
# Submit from the repository root:
#   sbatch ectoblob/run_ectoblob.sbatch
#SBATCH --time=4:00:00
#SBATCH --job-name=ectoblob
#SBATCH --partition=secondary
#SBATCH --account=campusclusterusers
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --array=0-4
#SBATCH --output=logs/ectoblob_%A_%a_log.txt

set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

cd "$PROJECT"
mkdir -p logs

# The Nth dataset directory under data/gene_trees, in sorted order.
DATASETS=(data/gene_trees/*)
DATASET=${DATASETS[$SLURM_ARRAY_TASK_ID]}
echo "=== task $SLURM_ARRAY_TASK_ID: $DATASET ==="

"$PYTHON" ectoblob/run_ectoblob.py \
    --config ectoblob/ectoblob_config.toml \
    --input "$DATASET" \
    --output "$DATASET"
