#!/usr/bin/env bash
# Copy the ASTRAL trees written next to the gene trees into a separate
# directory, keeping the subdirectory structure of the gene trees:
#
#   data/gene_trees/<dataset>/<rep>/astral.nwk
#     -> data/astral_trees/<dataset>/<rep>/astral.nwk
#
# Run after astral/run_all_astral.sh has finished (or at any point, to take a
# snapshot of the trees that exist so far):
#
#   bash astral/collect_astral_trees.sh
#   bash astral/collect_astral_trees.sh data/gene_trees data/astral_trees
set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
GENE_TREES=${1:-data/gene_trees}
OUT_ROOT=${2:-data/astral_trees}
OUT_NAME=${OUT_NAME:-astral.nwk}

cd "$PROJECT"
mkdir -p "$OUT_ROOT"
# Absolute paths, so the cd below cannot break the destination.
GENE_TREES=$(cd "$GENE_TREES" && pwd)
OUT_ROOT=$(cd "$OUT_ROOT" && pwd)

n=0
while IFS= read -r -d '' tree; do
    # find prints paths relative to $GENE_TREES; strip the leading "./"
    rel=${tree#./}
    dest=$OUT_ROOT/$rel
    mkdir -p "$(dirname "$dest")"
    cp -p "$GENE_TREES/$rel" "$dest"
    n=$((n + 1))
done < <(cd "$GENE_TREES" && find . -name "$OUT_NAME" -size +0 -print0)

echo "copied $n $OUT_NAME file(s) from $GENE_TREES to $OUT_ROOT"
