#!/usr/bin/env bash
# Infer a tree of blobs with ECToBlob for ONE replicate directory, passed
# as the first argument. Submitted once per replicate by
# ectoblob/run_all_ectoblob.sh, or directly:
#
#   sbatch ectoblob/ectoblob_one.sbatch data/gene_trees/miss_0_len_100/00
#SBATCH --time=04:00:00
#SBATCH --job-name=ectoblob-one
#SBATCH --partition=secondary
#SBATCH --account=campusclusterusers
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --output=logs/ectoblob_one_%j_log.txt

set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

if [ "$#" -lt 1 ]; then
    echo "usage: sbatch ectoblob/ectoblob_one.sbatch <replicate dir>" >&2
    exit 1
fi

cd "$PROJECT"
mkdir -p logs

echo "=== ECToBlob for $1 ==="
"$PYTHON" ectoblob/run_ectoblob.py \
    --config ectoblob/ectoblob_config.toml \
    --input "$1" \
    --output "$1"
