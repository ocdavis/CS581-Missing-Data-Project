#!/usr/bin/env bash
# Finish gene-tree inference for the subset replicates that are still
# missing their iqtree_500.nwk, as a SLURM job array: one task per missing
# replicate. The array is sized automatically at submit time to the
# replicates that are missing then (via redo_missing.py --list); each task
# runs infer_gene_trees.py on its replicate. Finished genes are skipped, so
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
#SBATCH --array=0-0%10
#SBATCH --output=logs/gene_trees_redo_%A_%a_log.txt

set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

cd "$PROJECT"
mkdir -p logs

# Size the array to the replicates that are missing right now: write the
# mapfile first, then rewrite this copy of the script with the real range
# and submit that. The submitted copy lives in logs/ and is not tracked.
MAPFILE=logs/redo_missing_mapfile.txt
"$PYTHON" dataset_maker/redo_missing.py --input data --list > "$MAPFILE"
N=$(wc -l < "$MAPFILE")
if [ "$N" -eq 0 ]; then
    echo "nothing to do: every replicate has its gene trees"
    exit 0
fi

SUBMIT=logs/redo_missing_array_${N}.sbatch
sed "s/^#SBATCH --array=0-0%10/#SBATCH --array=0-$((N - 1))%10/" \
    dataset_maker/redo_missing.sbatch > "$SUBMIT"
echo "submitting array of $N task(s): $SUBMIT"
sbatch "$SUBMIT"
exit 0

# --- the submitted array script continues below ----------------------------
# The array size was fixed at submit time; a task whose replicate was
# finished by another run in the meantime simply skips.

# The missing replicates, in the same order the array indices were sized
# against at submit time.
MAPFILE=logs/redo_missing_mapfile.txt
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
