#!/usr/bin/env bash
# Submit one ECToBlob job per replicate of the collected gene trees. Run
# this on the login node (not via sbatch): it lists the replicates without
# an iqtree_500.ectoblob.tre and sbatches ectoblob_one.sbatch for each, so
# every job is queued immediately and runs as resources allow.
#
# Usage (from the repository root):
#   bash ectoblob/run_all_ectoblob.sh                    # scan data/gene_trees
#   bash ectoblob/run_all_ectoblob.sh data/gene_trees/miss_0_len_100
set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

INPUT=${1:-data/gene_trees}
NAME=iqtree_500.ectoblob.tre

cd "$PROJECT"
mkdir -p logs

# Replicate directories (…/<dataset>/<rep>) without a ToB yet.
MAPFILE=logs/ectoblob_todo.txt
: > "$MAPFILE"
for rep in "$INPUT"/*/*/; do
    [ -f "$rep$NAME" ] || echo "$rep" >> "$MAPFILE"
done
N=$(wc -l < "$MAPFILE")
if [ "$N" -eq 0 ]; then
    echo "nothing to do: every replicate has its ToB"
    exit 0
fi

echo "submitting $N job(s) for the replicates listed in $MAPFILE"
i=0
while IFS= read -r rep; do
    i=$((i + 1))
    echo "[$i/$N] $rep"
    sbatch ectoblob/ectoblob_one.sbatch "$rep"
done < "$MAPFILE"
