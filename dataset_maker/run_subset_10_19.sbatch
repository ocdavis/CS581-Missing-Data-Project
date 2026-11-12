#!/usr/bin/env bash
# Generate the nine trimmed/missing-taxa datasets and infer IQ-TREE gene
# trees for replicates 10-19 of the five datasets of interest. All logic
# lives in dataset_maker/run_subset.py; this wrapper only sets the SLURM
# resources. See that script's docstring for the three steps (trim_phy,
# subset assembly, gene-tree inference).
#
# Submit from the repository root:
#   sbatch dataset_maker/run_subset_10_19.sbatch
#SBATCH --time=2-00:00:00
#SBATCH --job-name=gene-trees-10-19
#SBATCH --partition=cs
#SBATCH --account=cs-warnow
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g
#SBATCH --output=logs/gene_trees_10_19_%j_log.txt

set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

cd "$PROJECT"
mkdir -p logs

"$PYTHON" dataset_maker/run_subset.py 10-19 --threads "$SLURM_CPUS_PER_TASK"
