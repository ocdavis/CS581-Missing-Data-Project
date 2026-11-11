#!/usr/bin/env bash
# Time a full gene-tree inference on one replicate, to estimate how long
# the same calculation takes on every replicate.
#
# Replicate 00 of data/trimmed/len_500 is copied to scratch and given to
# infer_gene_trees.py (all 1000 genes, MFP selected on the first gene and
# reused for the rest, $SLURM_CPUS_PER_TASK concurrent IQ-TREE processes).
# The elapsed time is reported, along with the projection for all 150
# replicates (three lengths x 50) at this core count.
#
# Submit from the repository root:
#   sbatch dataset_maker/bench_mfp.sbatch
#SBATCH --time=04:00:00
#SBATCH --job-name=bench-replicate
#SBATCH --partition=secondary
#SBATCH --account=cs-warnow
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g
#SBATCH --output=logs/bench_mfp_%j_log.txt

set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python
REPLICATE=data/trimmed/len_500/00

cd "$PROJECT"
mkdir -p logs

SCRATCH=/scratch/quadros3/bench_replicate_$SLURM_JOB_ID/00
mkdir -p "$SCRATCH"
cp -r "$REPLICATE/seqs-500" "$SCRATCH/seqs-500"

echo "=== timing gene-tree inference on $REPLICATE (copy in $SCRATCH) ==="
T0=$(date +%s)
"$PYTHON" dataset_maker/infer_gene_trees.py \
    --model MFP
    --input "$SCRATCH" \
    --threads "$SLURM_CPUS_PER_TASK"
T1=$(date +%s)
ELAPSED=$((T1 - T0))

echo ""
echo "=== summary ==="
echo "$ELAPSED $SLURM_CPUS_PER_TASK" | awk '{
    printf "one replicate (%d genes, %d cores): %d min %d s\n", 1000, $2, $1 / 60, $1 % 60;
    printf "projected for 150 replicates:      %.1f h\n", 150 * $1 / 3600;
}'

cd /
rm -rf "$(dirname "$SCRATCH")"
echo "=== done ==="
