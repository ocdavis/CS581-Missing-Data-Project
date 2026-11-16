#!/usr/bin/env bash
# Estimate a TREE-QMC species tree for every gene-tree file under a
# directory. This is the base ("refinement") tree that both TOB-QMC and
# ECToBlob start from, so the two methods contract the same tree.
#
# For each <dir>/iqtree_500.nwk it runs
#   tree-qmc --root OUT -i <dir>/iqtree_500.nwk -o <dir>/tree-qmc.nwk
# and saves the messages in <dir>/tree-qmc.log. It then adds branch lengths
# to that tree, without changing it,
#   tree-qmc --root OUT --supportonly <dir>/tree-qmc.nwk \
#       -i <dir>/iqtree_500.nwk -o <dir>/tree-qmc.support.nwk
# and writes a plain-Newick copy as <dir>/tree-qmc.brlen.nwk:
#   tree-qmc.nwk          topology only
#   tree-qmc.support.nwk  + quartet support labels and branch lengths
#   tree-qmc.brlen.nwk    + branch lengths only (labels removed)
# Lengths are in coalescent units and exist for internal branches only (0
# where the gene trees do not favour the branch). Existing files are skipped.
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
# For another dataset, give its directory and, if they differ, the gene-tree
# file name, output name and outgroup:
#   GENES_NAME=genes.tre OUT_NAME=species.tre OUTGROUP=Homo_sapiens \
#       sbatch tob_qmc/run_base_tree.sbatch path/to/dataset
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
# Defaults; set these environment variables before submitting to change them
# (OUTGROUP= with no value leaves the tree unrooted).
export GENES_NAME=${GENES_NAME:-iqtree_500.nwk}
export OUT_NAME=${OUT_NAME:-tree-qmc.nwk}
export OUTGROUP=${OUTGROUP-OUT}

module load gcc/13.3.0

set -u
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PATH="${TREE_QMC_BUILD:-$PWD/../TREE-QMC/build}:$PATH"
command -v tree-qmc > /dev/null || { echo "tree-qmc not found on PATH"; exit 1; }

run_one() {
    genes=$1
    out=$(dirname "$genes")/$OUT_NAME
    stem=${out%.nwk}
    if [ -s "$out" ]; then
        echo "skip   $out"
    elif tree-qmc ${OUTGROUP:+--root "$OUTGROUP"} -i "$genes" -o "$out" --override \
            > "$stem.log" 2>&1 && [ -s "$out" ]; then
        echo "wrote  $out"
    else
        rm -f "$out"
        echo "FAILED $genes (see $stem.log)"
        return
    fi
    if [ -s "$stem.brlen.nwk" ]; then
        echo "skip   $stem.brlen.nwk"
    elif tree-qmc ${OUTGROUP:+--root "$OUTGROUP"} --supportonly "$out" -i "$genes" \
            -o "$stem.support.nwk" --override > "$stem.support.log" 2>&1 \
            && [ -s "$stem.support.nwk" ]; then
        # drop the quoted '[f1=...;EN=...]' labels, keep the :length values
        sed "s/'\[[^]]*\]'//g" "$stem.support.nwk" > "$stem.brlen.nwk"
        echo "wrote  $stem.brlen.nwk"
    else
        rm -f "$stem.support.nwk" "$stem.brlen.nwk"
        echo "FAILED branch lengths for $out (see $stem.support.log)"
    fi
}
export -f run_one

echo "=== base trees for $INPUT, ${SLURM_CPUS_PER_TASK:-1} at a time, $(command -v tree-qmc) ==="
find "$INPUT" -name "$GENES_NAME" -print0 | sort -z |
    xargs -0 -n 1 -P "${SLURM_CPUS_PER_TASK:-1}" bash -c 'run_one "$0"'

n_in=$(find "$INPUT" -name "$GENES_NAME" | wc -l)
n_out=$(find "$INPUT" -name "$OUT_NAME" -size +0 | wc -l)
n_len=$(find "$INPUT" -name "${OUT_NAME%.nwk}.brlen.nwk" -size +0 | wc -l)
echo "=== done: $n_out of $n_in base trees present, $n_len with branch lengths ==="
[ "$n_out" -eq "$n_in" ] && [ "$n_len" -eq "$n_in" ]
