#!/usr/bin/env bash
# Finish gene-tree inference for the subset replicates that are still
# missing their iqtree_500.nwk, as a SLURM job array: one task per missing
# replicate. dataset_maker/redo_missing.py --list finds the missing
# replicates (the array is sized to them at submit time); each task runs
# infer_gene_trees.py on its replicate. Finished genes are skipped, so
# partial replicates are completed rather than redone, and a task whose
# replicate was finished by another run exits immediately.
#
# Submit from the repository root:
#   sbatch dataset_maker/redo_missing.sbatch
#SBATCH --time=04:00:00
#SBATCH --job-name=gene-trees-redo
#SBATCH --partition=secondary
#SBATCH --account=campusclusterusers
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g
#SBATCH --array=0-39%10
#SBATCH --output=logs/gene_trees_redo_%A_%a_log.txt

set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

cd "$PROJECT"
mkdir -p logs

# The missing replicates, in the same order the array indices were sized
# against at submit time.
MAPFILE=logs/redo_missing_${SLURM_ARRAY_JOB_ID}.txt
"$PYTHON" dataset_maker/redo_missing.py --input data --list > "$MAPFILE"
N=$(wc -l < "$MAPFILE")
if [ "$SLURM_ARRAY_TASK_COUNT" -ne "$N" ]; then
    echo "array has $SLURM_ARRAY_TASK_COUNT tasks but $N replicates are missing" >&2
    echo "(resubmit to resize the array)" >&2
    exit 1
fi

REP=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MAPFILE")
echo "=== task $SLURM_ARRAY_TASK_ID: $REP ==="

# Skip if another task (or a previous run) already finished this replicate.
if [ -f "$REP/iqtree_500.nwk" ]; then
    echo "already done, skipping"
    exit 0
fi

"$PYTHON" dataset_maker/infer_gene_trees.py \
    --input "$REP" \
    --threads "$SLURM_CPUS_PER_TASK"
