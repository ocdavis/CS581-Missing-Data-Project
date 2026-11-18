#!/usr/bin/env bash
# Submit one gene-tree inference job per replicate for the extended
# conditions. Run this on the login node (not via sbatch): it lists the
# replicates that still lack their iqtree_500.nwk and sbatches
# gene_trees_one.sbatch for each, so every job is queued immediately and runs
# as resources allow.
#
# Nine datasets x 50 replicates = 450 jobs, written to the usual layout:
#
#   data/gene_trees/<dataset>/{00..49}/iqtree_500.nwk
#
# covering
#   miss_0_len_50_files_100, miss_0_len_25_files_100      (shorter sequences)
#   miss_50_len_500_files_95, miss_50_len_500_files_99    (50% taxa dropped)
#   miss_80_len_500_files_{100,50,80,95,99}               (80% taxa dropped)
#
# Settings match the existing gene trees: model GTR+G, seed 581, 32
# concurrent IQ-TREE processes per replicate, work folders removed on
# success.
#
# The account caps *running* jobs at 40, so only 40 run at a time and the
# rest wait in the queue; submission itself is not limited. A run of this
# size therefore takes roughly 450/40 x 10 min, and re-running the script
# later picks up anything that failed or has not started yet. Replicates
# whose iqtree_500.nwk already exists are skipped, so relaunching is safe.
#
# Usage (from the repository root):
#   bash dataset_maker/run_extended_trees.sh
#   bash dataset_maker/run_extended_trees.sh miss_0_len_25_files_100   # just one
set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
DATASETS=${*:-"miss_0_len_50_files_100 miss_0_len_25_files_100 \
miss_50_len_500_files_95 miss_50_len_500_files_99 \
miss_80_len_500_files_100 miss_80_len_500_files_50 miss_80_len_500_files_80 \
miss_80_len_500_files_95 miss_80_len_500_files_99"}
NAME=iqtree_500.nwk
TRIMMED=data/trimmed
GENE_TREES=data/gene_trees

cd "$PROJECT"
mkdir -p logs

# Absolute replicate paths, one per line.
MAPFILE=logs/extended_trees_todo.txt
: > "$MAPFILE"
for ds in $DATASETS; do
    if [ ! -d "$TRIMMED/$ds" ]; then
        echo "error: $TRIMMED/$ds not found;" \
             "run 'sbatch dataset_maker/run_extended_datasets.sbatch' first" >&2
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
echo "(the account allows 40 running at once, so most will wait in the queue)"
i=0
failed=0
while IFS= read -r rep; do
    i=$((i + 1))
    ds=$(basename "$(dirname "$rep")")
    out="$GENE_TREES/$ds/$(basename "$rep")"
    echo "[$i/$N] $rep -> $out"
    # Keep going if a single submission is rejected, and report it at the end
    # instead of aborting the batch partway through.
    if ! sbatch dataset_maker/gene_trees_one.sbatch "$rep" "$out"; then
        echo "  WARNING: submission failed for $rep" >&2
        failed=$((failed + 1))
    fi
done < "$MAPFILE"

if [ "$failed" -gt 0 ]; then
    echo "$failed of $N submission(s) failed; re-run this script to retry" >&2
    exit 1
fi
echo "submitted $N job(s)"
