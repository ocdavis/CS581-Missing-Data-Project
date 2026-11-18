#!/usr/bin/env bash
# Submit one ECToBlob job per replicate of the partial-missingness datasets.
# Run this on the login node (not via sbatch): it lists the replicates that
# still lack their iqtree_500.ectoblob.tre and sbatches ectoblob_one.sbatch
# for each, so all 100 jobs are queued immediately and run as resources
# allow.
#
# Each job reads, from the replicate directory:
#
#   iqtree_500.nwk        the inferred gene trees (one per line)
#   tree-qmc.brlen.nwk    the TREE-QMC species tree used as the starting
#                         tree / resolution of the true tree of blobs
#
# and writes next to them, matching the layout of the other datasets:
#
#   data/gene_trees/miss_50_len_500_files_50/{00..49}/iqtree_500.ectoblob.tre
#   data/gene_trees/miss_50_len_500_files_80/{00..49}/iqtree_500.ectoblob.tre
#
# plus a matching .ectoblob.log with R's messages.
#
# Replicates that already have a ToB are skipped, so an interrupted batch can
# simply be relaunched and it will also retry any replicate whose job failed.
# Replicates whose gene trees or starting tree are not there yet are reported
# and skipped rather than submitted to fail, so this is safe to run before
# the TREE-QMC trees exist and again once they do.
#
# Usage (from the repository root):
#   bash ectoblob/run_partial_miss_ectoblob.sh
#   bash ectoblob/run_partial_miss_ectoblob.sh miss_50_len_500_files_50
set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
DATASETS=${*:-"miss_50_len_500_files_50 miss_50_len_500_files_80"}
GENE_TREES=data/gene_trees
GENE_TREE_NAME=iqtree_500.nwk
START_TREE_NAME=tree-qmc.brlen.nwk
TOB_NAME=iqtree_500.ectoblob.tre

cd "$PROJECT"
mkdir -p logs

# Absolute replicate paths, one per line.
MAPFILE=logs/partial_miss_ectoblob_todo.txt
: > "$MAPFILE"
n_done=0
n_no_genes=0
n_no_start=0
for ds in $DATASETS; do
    if [ ! -d "$GENE_TREES/$ds" ]; then
        echo "error: $GENE_TREES/$ds not found;" \
             "run dataset_maker/run_partial_miss_trees.sh first" >&2
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
while IFS= read -r rep; do
    i=$((i + 1))
    echo "[$i/$N] $rep"
    sbatch ectoblob/ectoblob_one.sbatch "$rep"
done < "$MAPFILE"
