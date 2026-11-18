#!/usr/bin/env bash
# Infer gene trees for ONE replicate directory, passed as the first
# argument. Submitted once per replicate by redo_missing.sh and
# run_partial_miss_trees.sh.
#
# With one argument the trees are written next to the replicate's seqs-500
# folder. With a second argument they are written to that directory instead
# (as <output dir>/iqtree_500.nwk) and the IQ-TREE work folder is deleted
# afterwards, which keeps data/trimmed free of intermediates.
#
# Usage (from the repository root):
#   sbatch dataset_maker/gene_trees_one.sbatch data/subset_20_24/miss_0_len_500/24
#   sbatch dataset_maker/gene_trees_one.sbatch \
#       data/trimmed/miss_50_len_500_files_50/00 \
#       data/gene_trees/miss_50_len_500_files_50/00
#SBATCH --time=20:00
#SBATCH --job-name=gene-trees-one
#SBATCH --partition=scavenger,secondary,cs
#SBATCH --account=campusclusterusers
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g
#SBATCH --output=logs/gene_trees_one_%j_log.txt

set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

if [ "$#" -lt 1 ]; then
    echo "usage: sbatch dataset_maker/gene_trees_one.sbatch <replicate dir> [output dir]" >&2
    exit 1
fi

cd "$PROJECT"
mkdir -p logs

echo "=== iqtree: $(command -v iqtree || echo MISSING) ==="
echo "=== inferring gene trees for $1 ==="
if [ "$#" -ge 2 ]; then
    echo "=== writing $2/iqtree_500.nwk (work folder removed on success) ==="
    "$PYTHON" dataset_maker/infer_gene_trees.py \
        --input "$1" \
        --output "$2" \
        --threads "$SLURM_CPUS_PER_TASK" \
        --clean
else
    "$PYTHON" dataset_maker/infer_gene_trees.py \
        --input "$1" \
        --threads "$SLURM_CPUS_PER_TASK"
fi
