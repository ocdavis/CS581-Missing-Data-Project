#!/usr/bin/env bash
# Submit one ASTRAL job per dataset directory under data/gene_trees.
# Run this on the login node (not via sbatch): it lists the directories that
# still lack some astral.nwk and sbatches astral/astral_one.sbatch for each,
# so every job is queued immediately and runs as resources allow.
#
# Each job runs astral/run_astral.py on its directory, which finds every
# iqtree_500.nwk below it and writes astral.nwk (plus astral.log) next to it:
#
#   data/gene_trees/<dataset>/<rep>/iqtree_500.nwk
#   data/gene_trees/<dataset>/<rep>/astral.nwk
#
# To gather the trees in a separate directory afterwards, keeping the
# subdirectory structure, run astral/collect_astral_trees.sh.
#
# Directories that already have every tree are skipped, so an interrupted
# batch can simply be relaunched and it will also retry directories whose job
# failed.
#
# Usage (from the repository root):
#   bash astral/run_all_astral.sh
#   bash astral/run_all_astral.sh miss_0_len_500 miss_50_len_500
set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
GENE_TREES=data/gene_trees
GENES_NAME=${GENES_NAME:-iqtree_500.nwk}
OUT_NAME=${OUT_NAME:-astral.nwk}

cd "$PROJECT"
mkdir -p logs

if [ "$#" -gt 0 ]; then
    DATASETS=$*
else
    DATASETS=$(cd "$GENE_TREES" && ls -d */ | sed 's#/##')
fi

# Absolute dataset paths, one per line.
MAPFILE=logs/astral_todo.txt
: > "$MAPFILE"
n_done=0
n_no_genes=0
for ds in $DATASETS; do
    dir=$GENE_TREES/$ds
    if [ ! -d "$dir" ]; then
        echo "error: $dir not found" >&2
        exit 1
    fi
    n_in=$(find "$dir" -name "$GENES_NAME" | wc -l)
    n_out=$(find "$dir" -name "$OUT_NAME" -size +0 | wc -l)
    if [ "$n_in" -eq 0 ]; then
        echo "skip, no $GENES_NAME: $dir" >&2
        n_no_genes=$((n_no_genes + 1))
    elif [ "$n_out" -ge "$n_in" ]; then
        n_done=$((n_done + 1))
    else
        echo "$dir" >> "$MAPFILE"
    fi
done

N=$(wc -l < "$MAPFILE")
if [ "$N" -eq 0 ]; then
    echo "nothing to submit: $n_done dataset(s) already have all $OUT_NAME," \
         "$n_no_genes missing $GENES_NAME"
    exit 0
fi

echo "submitting $N job(s) for the directories listed in $MAPFILE"
if [ "$n_done" -gt 0 ]; then
    echo "  ($n_done already done)"
fi
if [ "$n_no_genes" -gt 0 ]; then
    echo "  ($n_no_genes skipped: no $GENES_NAME)"
fi

i=0
failed=0
while IFS= read -r dir; do
    i=$((i + 1))
    echo "[$i/$N] $dir"
    # Keep going if a single submission is rejected, and report it at the end
    # instead of aborting the batch partway through.
    if ! sbatch astral/astral_one.sbatch "$dir"; then
        echo "  WARNING: submission failed for $dir" >&2
        failed=$((failed + 1))
    fi
done < "$MAPFILE"

if [ "$failed" -gt 0 ]; then
    echo "$failed of $N submission(s) failed; re-run this script to retry" >&2
    exit 1
fi
echo "submitted $N job(s)"
