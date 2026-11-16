#!/usr/bin/env bash
# Estimate a TREE-QMC species tree for every gene-tree file under a
# directory. This is the base ("refinement") tree that both TOB-QMC and
# ECToBlob start from, so the two methods contract the same tree.
#
# For each <dir>/iqtree_500.nwk it runs
#   tree-qmc --root OUT -i <dir>/iqtree_500.nwk -o <dir>/tree-qmc.nwk
# and saves the messages in <dir>/tree-qmc.log. Existing trees are skipped.
#
# Expects tree-qmc in a TREE-QMC checkout next to this repository
# (../TREE-QMC/build/tree-qmc); set TREE_QMC_BUILD to use another location.
#
# Submit from the repository root (input directory defaults to
# data/gene_trees):
#   mkdir -p logs
#   sbatch tob_qmc/run_base_tree.sbatch data/gene_trees
# Options given before the script name override the #SBATCH lines, e.g. a
# different account, or either partition with a time limit that fits both:
#   sbatch --account=NAME --partition=cs,secondary --time=01:00:00 \
#       tob_qmc/run_base_tree.sbatch data/gene_trees
# or run without SLURM on a small directory:
#   bash tob_qmc/run_base_tree.sbatch data/gene_trees/miss_0_len_500/00
#SBATCH --time=1-00:00:00
#SBATCH --job-name=base-trees
#SBATCH --partition=cs
#SBATCH --account=cs-melkebir
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --output=logs/base_trees_%j_log.txt

INPUT=${1:-data/gene_trees}
export GENES_NAME=iqtree_500.nwk OUT_NAME=tree-qmc.nwk OUTGROUP=OUT

module load gcc/13.3.0

set -u
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PATH="${TREE_QMC_BUILD:-$PWD/../TREE-QMC/build}:$PATH"
command -v tree-qmc > /dev/null || { echo "tree-qmc not found on PATH"; exit 1; }

run_one() {
    genes=$1
    out=$(dirname "$genes")/$OUT_NAME
    if [ -s "$out" ]; then
        echo "skip   $out"
    elif tree-qmc --root "$OUTGROUP" -i "$genes" -o "$out" --override \
            > "${out%.nwk}.log" 2>&1 && [ -s "$out" ]; then
        echo "wrote  $out"
    else
        rm -f "$out"
        echo "FAILED $genes (see ${out%.nwk}.log)"
    fi
}
export -f run_one

echo "=== base trees for $INPUT, ${SLURM_CPUS_PER_TASK:-1} at a time, $(command -v tree-qmc) ==="
find "$INPUT" -name "$GENES_NAME" -print0 | sort -z |
    xargs -0 -n 1 -P "${SLURM_CPUS_PER_TASK:-1}" bash -c 'run_one "$0"'

n_in=$(find "$INPUT" -name "$GENES_NAME" | wc -l)
n_out=$(find "$INPUT" -name "$OUT_NAME" -size +0 | wc -l)
echo "=== done: $n_out of $n_in base trees present ==="
[ "$n_out" -eq "$n_in" ]
