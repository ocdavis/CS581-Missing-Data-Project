#!/usr/bin/env bash
# Submit one gene-tree inference job per replicate of the partial-missingness
# datasets. Run this on the login node (not via sbatch): it lists the
# replicates that still lack their iqtree_500.nwk and sbatches
# gene_trees_one.sbatch for each, so all 100 jobs are queued immediately and
# run as resources allow.
#
# Result, matching the layout of the other datasets:
#
#   data/gene_trees/miss_50_files_50_len_500/{00..49}/iqtree_500.nwk
#   data/gene_trees/miss_50_files_80_len_500/{00..49}/iqtree_500.nwk
#
# Each file holds one ML tree per line, in gene order, for the 1000 genes of
# that replicate.
#
# Replicates that already have a gene-tree file are skipped, so an
# interrupted batch can simply be relaunched (re-running also picks up any
# replicate whose job failed).
#
# Usage (from the repository root):
#   bash dataset_maker/run_partial_miss_trees.sh
#   bash dataset_maker/run_partial_miss_trees.sh miss_50_files_50_len_500
set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
DATASETS=${*:-"miss_50_files_50_len_500 miss_50_files_80_len_500"}
NAME=iqtree_500.nwk
TRIMMED=data/trimmed
GENE_TREES=data/gene_trees

cd "$PROJECT"
mkdir -p logs

# Absolute replicate paths, one per line.
MAPFILE=logs/partial_miss_trees_todo.txt
: > "$MAPFILE"
for ds in $DATASETS; do
    if [ ! -d "$TRIMMED/$ds" ]; then
        echo "error: $TRIMMED/$ds not found; run run_partial_miss.sbatch first" >&2
        exit 1
    fi
    for rep in "$TRIMMED/$ds"/*/; do
        rep=${rep%/}
        [ -f "$GENE_TREES/$ds/$(basename "$rep")/$NAME" ] || echo "$rep" >> "$MAPFILE"
    done
done

N=$(wc -l < "$MAPFILE")
if [ "$N" -eq 0 ]; then
    echo "nothing to do: every replicate has its gene trees"
    exit 0
fi

echo "submitting $N job(s) for the replicates listed in $MAPFILE"
i=0
while IFS= read -r rep; do
    i=$((i + 1))
    ds=$(basename "$(dirname "$rep")")
    out="$GENE_TREES/$ds/$(basename "$rep")"
    echo "[$i/$N] $rep -> $out"
    sbatch dataset_maker/gene_trees_one.sbatch "$rep" "$out"
done < "$MAPFILE"
