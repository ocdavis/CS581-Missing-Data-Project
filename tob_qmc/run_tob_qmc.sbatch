#!/usr/bin/env bash
# Estimate a tree of blobs with TOB-QMC for every gene-tree file under a
# directory, starting from the TREE-QMC base tree made by
# run_base_tree.sbatch (the same starting tree ECToBlob uses).
#
# For each folder holding iqtree_500.nwk and tree-qmc.nwk it runs
#   1. p-value search around every branch of the base tree (slow):
#        tree-qmc --blobsearchonly tree-qmc.nwk --store_pvalue \
#            --iter_limit_blob 0 -i iqtree_500.nwk -o tob-qmc.pvalues.nwk
#   2. contraction at the chosen significance levels (instant):
#        tree-qmc --alpha 1e-7 --beta 0.95 --load_pvalue \
#            -i tob-qmc.pvalues.nwk -o tob-qmc.nwk
# with the messages in tob-qmc.search.log and tob-qmc.contract.log. Finished
# steps are skipped, so an interrupted job can be resubmitted, and deleting
# only tob-qmc.nwk redoes the contraction (e.g. with another ALPHA/BETA)
# without repeating the search.
#
# Needs tree-qmc built WITH tree-of-blobs support (cmake + R, see README), in
# a TREE-QMC-tob checkout next to this repository
# (../TREE-QMC-tob/build/tree-qmc); set TREE_QMC_BUILD for another location.
#
# Submit from the repository root (input directory defaults to
# data/gene_trees):
#   mkdir -p logs
#   sbatch tob_qmc/run_tob_qmc.sbatch data/gene_trees
# Options before the script name override the #SBATCH lines, e.g.
#   sbatch --account=NAME --partition=cs,secondary --time=04:00:00 \
#       tob_qmc/run_tob_qmc.sbatch data/gene_trees
# or run without SLURM on one folder:
#   bash tob_qmc/run_tob_qmc.sbatch data/gene_trees/miss_0_len_500/00
#SBATCH --time=2-00:00:00
#SBATCH --job-name=tob-qmc
#SBATCH --partition=cs
#SBATCH --account=cs-melkebir
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --output=logs/tob_qmc_%j_log.txt

INPUT=${1:-data/gene_trees}
# Defaults; set these environment variables before submitting to change them.
export GENES_NAME=${GENES_NAME:-iqtree_500.nwk}
export BASE_NAME=${BASE_NAME:-tree-qmc.nwk}
export PVAL_NAME=${PVAL_NAME:-tob-qmc.pvalues.nwk}
export TOB_NAME=${TOB_NAME:-tob-qmc.nwk}
export ALPHA=${ALPHA:-1e-7}
export BETA=${BETA:-0.95}
# Search options: exhaustive by default; SEARCH= (empty) for TOB-QMC's
# default heuristic, SEARCH=--3f1a for 3-fix-1-alter.
export SEARCH=${SEARCH---iter_limit_blob 0}

set -u
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export PATH="${TREE_QMC_BUILD:-$PWD/../TREE-QMC-tob/build}:$PATH"
command -v tree-qmc > /dev/null || { echo "tree-qmc not found on PATH"; exit 1; }

run_one() {
    dir=$(dirname "$1")
    genes=$1 base=$dir/$BASE_NAME pval=$dir/$PVAL_NAME tob=$dir/$TOB_NAME
    stem=${tob%.nwk}
    if [ -s "$tob" ]; then
        echo "skip   $tob"
        return
    fi
    if [ ! -s "$base" ]; then
        echo "FAILED $dir: no base tree $BASE_NAME (run run_base_tree.sbatch first)"
        return
    fi
    if [ ! -s "$pval" ]; then
        start=$SECONDS
        # shellcheck disable=SC2086
        if ! tree-qmc --blobsearchonly "$base" --store_pvalue $SEARCH \
                -i "$genes" -o "$pval.tmp" --override > "$stem.search.log" 2>&1 \
                || [ ! -s "$pval.tmp" ]; then
            rm -f "$pval.tmp"
            echo "FAILED $dir: p-value search (see $stem.search.log)"
            return
        fi
        mv "$pval.tmp" "$pval"
        echo "search $pval ($((SECONDS - start)) s)"
    fi
    if tree-qmc --alpha "$ALPHA" --beta "$BETA" --load_pvalue \
            -i "$pval" -o "$tob" --override > "$stem.contract.log" 2>&1 \
            && [ -s "$tob" ]; then
        echo "wrote  $tob"
    else
        rm -f "$tob"
        echo "FAILED $dir: contraction (see $stem.contract.log)"
    fi
}
export -f run_one

echo "=== TOB-QMC for $INPUT: alpha $ALPHA, beta $BETA, search '${SEARCH:-default}', ${SLURM_CPUS_PER_TASK:-1} at a time, $(command -v tree-qmc) ==="
find "$INPUT" -name "$GENES_NAME" -print0 | sort -z |
    xargs -0 -n 1 -P "${SLURM_CPUS_PER_TASK:-1}" bash -c 'run_one "$0"'

n_in=$(find "$INPUT" -name "$GENES_NAME" | wc -l)
n_out=$(find "$INPUT" -name "$TOB_NAME" -size +0 | wc -l)
echo "=== done: $n_out of $n_in trees of blobs present ==="
[ "$n_out" -eq "$n_in" ]
