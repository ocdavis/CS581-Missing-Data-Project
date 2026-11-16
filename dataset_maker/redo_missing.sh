#!/usr/bin/env bash
# Submit one gene-tree inference job per replicate that is still missing
# its iqtree_500.nwk. Run this on the login node (not via sbatch): it
# lists the missing replicates and sbatches gene_trees_one.sbatch for
# each, so every job is queued immediately and runs as resources allow.
#
# Usage (from the repository root):
#   bash dataset_maker/redo_missing.sh                 # scan data/
#   bash dataset_maker/redo_missing.sh data/subset_20_24
set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

INPUT=${1:-data}

cd "$PROJECT"
mkdir -p logs

# Absolute replicate paths, one per line.
MAPFILE=logs/redo_missing_mapfile.txt
"$PYTHON" dataset_maker/redo_missing.py --input "$INPUT" --list > "$MAPFILE"
N=$(wc -l < "$MAPFILE")
if [ "$N" -eq 0 ]; then
    echo "nothing to do: every replicate has its gene trees"
    exit 0
fi

echo "submitting $N job(s) for the replicates listed in $MAPFILE"
i=0
while IFS= read -r rep; do
    i=$((i + 1))
    echo "[$i/$N] $rep"
    sbatch dataset_maker/gene_trees_one.sbatch "$rep"
done < "$MAPFILE"
