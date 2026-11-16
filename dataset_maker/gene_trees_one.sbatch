#!/usr/bin/env bash
# Infer gene trees for ONE replicate directory, passed as the first
# argument. Submitted once per missing replicate by redo_missing.sh.
#
# Usage (from the repository root):
#   sbatch dataset_maker/gene_trees_one.sbatch data/subset_20_24/miss_0_len_500/24
#SBATCH --time=20:00
#SBATCH --job-name=gene-trees-one
#SBATCH --partition=scavenger
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
    echo "usage: sbatch dataset_maker/gene_trees_one.sbatch <replicate dir>" >&2
    exit 1
fi

cd "$PROJECT"
mkdir -p logs

echo "=== inferring gene trees for $1 ==="
"$PYTHON" dataset_maker/infer_gene_trees.py \
    --input "$1" \
    --threads "$SLURM_CPUS_PER_TASK"
