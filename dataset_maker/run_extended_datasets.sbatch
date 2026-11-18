#!/usr/bin/env bash
# Generate the extended set of gene-tree-estimation-error conditions: shorter
# sequences, and 80% missing taxa, at the same full-length layout as the
# existing datasets.
#
# Nine datasets are written under data/trimmed (all with seed 581, so the taxa
# dropped from a given gene match the other missing-taxa datasets):
#
#   shorter sequences, all taxa present
#     miss_0_len_50_files_100               50 sites
#     miss_0_len_25_files_100               25 sites
#
#   50% of taxa dropped, full length
#     miss_50_len_500_files_95              95% of genes affected
#     miss_50_len_500_files_99              99% of genes affected
#
#   80% of taxa dropped, full length
#     miss_80_len_500_files_100             every gene affected
#     miss_80_len_500_files_50              50% of genes affected
#     miss_80_len_500_files_80              80% of genes affected
#     miss_80_len_500_files_95              95% of genes affected
#     miss_80_len_500_files_99              99% of genes affected
#
# With 50 ingroup taxa, 50% removes 25 (26 rows left) and 80% removes 40
# (11 rows left); the outgroup is never removed. Existing datasets are
# untouched: trim_phy.py only writes the directories named above.
#
# The three groups cannot be expressed in one trim_phy config file, because a
# config crosses every length with every missing and file percentage; each
# group is therefore a separate invocation. No config file is used at all —
# the parameters are on the command lines below.
#
# Submit from the repository root:
#   sbatch dataset_maker/run_extended_datasets.sbatch
#SBATCH --time=04:00:00
#SBATCH --job-name=trim-extended
#SBATCH --partition=secondary
#SBATCH --account=campusclusterusers
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g
#SBATCH --output=logs/trim_extended_%j_log.txt

set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

cd "$PROJECT"
mkdir -p logs

INPUT=data/camus-dataset/n50
OUTPUT=data/trimmed
SEED=581
W="$SLURM_CPUS_PER_TASK"

echo "=== 1/3: shortened sequences, all taxa ==="
"$PYTHON" dataset_maker/trim_phy.py \
    --input "$INPUT" --output "$OUTPUT" \
    --lengths 50 25 --missing 0 \
    --seed "$SEED" --workers "$W"

echo "=== 2/3: 50% of taxa dropped from 95% / 99% of genes ==="
"$PYTHON" dataset_maker/trim_phy.py \
    --input "$INPUT" --output "$OUTPUT" \
    --lengths 500 --missing 50 --file-percents 95 99 \
    --seed "$SEED" --workers "$W"

echo "=== 3/3: 80% of taxa dropped from 50/80/95/99/100% of genes ==="
"$PYTHON" dataset_maker/trim_phy.py \
    --input "$INPUT" --output "$OUTPUT" \
    --lengths 500 --missing 80 --file-percents 50 80 95 99 100 \
    --seed "$SEED" --workers "$W"

# Informational only: trim_phy already failed above if anything went wrong.
echo "=== datasets written under $OUTPUT ==="
ls -d "$OUTPUT"/miss_0_len_50_files_100 "$OUTPUT"/miss_0_len_25_files_100 \
      "$OUTPUT"/miss_50_len_500_files_95 "$OUTPUT"/miss_50_len_500_files_99 \
      "$OUTPUT"/miss_80_len_500_files_100 "$OUTPUT"/miss_80_len_500_files_50 \
      "$OUTPUT"/miss_80_len_500_files_80 "$OUTPUT"/miss_80_len_500_files_95 \
      "$OUTPUT"/miss_80_len_500_files_99 2>/dev/null || true
