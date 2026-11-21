#!/usr/bin/env bash
# Infer ASTRAL species trees for every gene-tree file under ONE directory,
# passed as the first argument. Submitted once per dataset directory by
# astral/run_all_astral.sh, or directly:
#
#   sbatch astral/astral_one.sbatch data/gene_trees/miss_0_len_500
#
# Extra arguments are passed to run_astral.py, e.g. to write the trees under
# data/astral_trees/<dataset>/ instead of next to the gene trees:
#
#   sbatch astral/astral_one.sbatch data/gene_trees/miss_0_len_500 \
#       --output data/astral_trees/miss_0_len_500
#SBATCH --time=4:00:00
#SBATCH --job-name=astral-one
#SBATCH --partition=secondary,cs,scavenger
#SBATCH --account=campusclusterusers
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --output=logs/astral_one_%j_log.txt

set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

if [ "$#" -lt 1 ]; then
    echo "usage: sbatch astral/astral_one.sbatch <directory> [run_astral.py options]" >&2
    exit 1
fi

cd "$PROJECT"
mkdir -p logs

echo "=== ASTRAL for $1 ==="
"$PYTHON" astral/run_astral.py "$@"
