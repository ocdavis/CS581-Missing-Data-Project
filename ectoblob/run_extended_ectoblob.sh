#!/usr/bin/env bash
# Submit one ECToBlob job per replicate of the extended gene-tree datasets.
# Run this on the login node (not via sbatch): it lists the replicates that
# still lack their iqtree_500.ectoblob.tre and sbatches ectoblob_one.sbatch
# for each, so every job is queued immediately and runs as resources allow.
#
# Each job reads, from the replicate directory:
#
#   iqtree_500.nwk        the inferred gene trees (one per line)
#   tree-qmc.brlen.nwk    the TREE-QMC species tree used as the starting
#                         tree / resolution of the true tree of blobs
#
# and writes next to them, matching the layout of the other datasets:
#
#   data/gene_trees/<dataset>/{00..49}/iqtree_500.ectoblob.tre
#
# plus a matching .ectoblob.log with R's messages.
#
# Eight datasets x 50 replicates = 400 jobs:
#
#   miss_0_len_50_files_100, miss_0_len_25_files_100      (shorter sequences)
#   miss_50_len_500_files_95, miss_50_len_500_files_99    (50% taxa dropped)
#   miss_80_len_500_files_{50,80,95,99}                   (80% taxa dropped)
#
# miss_80_len_500_files_100 is deliberately NOT included. Every gene in that
# dataset has 80% of its taxa dropped, and different genes drop different
# taxa, so no single gene tree contains the full taxon set. ECToBlob takes
# the taxon set from the gene trees and requires it to match the starting
# tree, so every replicate there would fail; there is nothing to infer.
# The other datasets all contain complete gene trees (1% to 100% of genes),
# which ectoblob.R promotes to the front so the taxon set matches.
#
# Replicates that already have a ToB are skipped, so an interrupted batch can
# simply be relaunched and it will also retry any replicate whose job failed.
# Replicates whose gene trees or starting tree are missing are reported and
# skipped rather than submitted to fail.
#
# Usage (from the repository root):
#   bash ectoblob/run_extended_ectoblob.sh
#   bash ectoblob/run_extended_ectoblob.sh miss_80_len_500_files_99
set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
DATASETS=${*:-"miss_0_len_50_files_100 miss_0_len_25_files_100 \
miss_50_len_500_files_95 miss_50_len_500_files_99 \
miss_80_len_500_files_50 miss_80_len_500_files_80 \
miss_80_len_500_files_95 miss_80_len_500_files_99"}
GENE_TREES=data/gene_trees
GENE_TREE_NAME=iqtree_500.nwk
START_TREE_NAME=tree-qmc.brlen.nwk
TOB_NAME=iqtree_500.ectoblob.tre

cd "$PROJECT"
mkdir -p logs

# Absolute replicate paths, one per line.
MAPFILE=logs/extended_ectoblob_todo.txt
: > "$MAPFILE"
n_done=0
n_no_genes=0
n_no_start=0
for ds in $DATASETS; do
    if [ ! -d "$GENE_TREES/$ds" ]; then
        echo "error: $GENE_TREES/$ds not found;" \
             "run dataset_maker/run_extended_trees.sh first" >&2
        exit 1
    fi
    for rep in "$GENE_TREES/$ds"/*/; do
        rep=${rep%/}
        if [ -f "$rep/$TOB_NAME" ]; then
            n_done=$((n_done + 1))
        elif [ ! -f "$rep/$GENE_TREE_NAME" ]; then
            echo "skip, no $GENE_TREE_NAME: $rep" >&2
            n_no_genes=$((n_no_genes + 1))
        elif [ ! -f "$rep/$START_TREE_NAME" ]; then
            echo "skip, no $START_TREE_NAME: $rep" >&2
            n_no_start=$((n_no_start + 1))
        else
            echo "$rep" >> "$MAPFILE"
        fi
    done
done

N=$(wc -l < "$MAPFILE")
if [ "$N" -eq 0 ]; then
    echo "nothing to submit: $n_done replicate(s) already have a ToB," \
         "$n_no_start missing $START_TREE_NAME," \
         "$n_no_genes missing $GENE_TREE_NAME"
    exit 0
fi

echo "submitting $N job(s) for the replicates listed in $MAPFILE"
echo "(the account allows 40 running at once, so most will wait in the queue)"
if [ "$n_done" -gt 0 ]; then
    echo "  ($n_done already done)"
fi
if [ "$n_no_start" -gt 0 ]; then
    echo "  ($n_no_start skipped: no $START_TREE_NAME yet)"
fi
if [ "$n_no_genes" -gt 0 ]; then
    echo "  ($n_no_genes skipped: no $GENE_TREE_NAME yet)"
fi

i=0
failed=0
while IFS= read -r rep; do
    i=$((i + 1))
    echo "[$i/$N] $rep"
    # Keep going if a single submission is rejected, and report it at the end
    # instead of aborting the batch partway through.
    if ! sbatch ectoblob/ectoblob_one.sbatch "$rep"; then
        echo "  WARNING: submission failed for $rep" >&2
        failed=$((failed + 1))
    fi
done < "$MAPFILE"

if [ "$failed" -gt 0 ]; then
    echo "$failed of $N submission(s) failed; re-run this script to retry" >&2
    exit 1
fi
echo "submitted $N job(s)"
