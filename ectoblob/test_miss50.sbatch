#!/usr/bin/env bash
# Test ECToBlob on a full miss_50 replicate with a complete first gene tree:
# replicate 00 of miss_50_len_500 is copied to scratch, its first gene tree
# is replaced with the first (complete, 51-taxon) gene tree of miss_0_len_500
# replicate 00, and ectoblob.R is run with the same parameters as the real
# pipeline (from ectoblob_config.toml). The point is to check whether the
# mul-path star-tree crash reproduces with the full 1000 gene trees.
#
# Submit from the repository root:
#   sbatch ectoblob/test_miss50.sbatch
#SBATCH --time=04:00:00
#SBATCH --job-name=ectoblob-test
#SBATCH --partition=secondary
#SBATCH --account=campusclusterusers
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --output=logs/ectoblob_test_%j_log.txt

set -eu

PROJECT=/projects/illinois/eng/cs/warnow/quadros3/classes/cs581/cs581-final-project
PYTHON=/u/quadros3/miniconda3/envs/cs581-final-project/bin/python

cd "$PROJECT"
mkdir -p logs

# --- assemble the test replicate in scratch -------------------------------
SCRATCH=/scratch/quadros3/ectoblob_miss50_test_$SLURM_JOB_ID
mkdir -p "$SCRATCH"

# The miss_50 replicate: gene trees + starting tree.
cp data/gene_trees/miss_50_len_500/00/iqtree_500.nwk "$SCRATCH/genes.nwk"
cp data/gene_trees/miss_50_len_500/00/tree-qmc.brlen.nwk "$SCRATCH/start.nwk"

# Replace the first gene tree with the first (complete) gene tree of the
# miss_0 replicate, so quartetTable builds its columns over all 51 taxa.
head -n 1 data/gene_trees/miss_0_len_500/00/iqtree_500.nwk > "$SCRATCH/genes.nwk.new"
tail -n +2 "$SCRATCH/genes.nwk" >> "$SCRATCH/genes.nwk.new"
mv "$SCRATCH/genes.nwk.new" "$SCRATCH/genes.nwk"

echo "=== test: miss_50_len_500/00 with complete first gene tree ==="
echo "gene trees: $(wc -l < "$SCRATCH/genes.nwk")"
echo "scratch: $SCRATCH"

# --- run ECToBlob with the real pipeline's parameters ---------------------
# Read the parameters from ectoblob_config.toml so this test always matches
# the production configuration.
eval "$("$PYTHON" - <<'EOF'
import tomllib
cfg = tomllib.load(open("ectoblob/ectoblob_config.toml", "rb"))
print(f"ALPHA={cfg.get('alpha', 0.05):g}")
print(f"BETA={cfg.get('beta', 0.8):g}")
print(f"QTYPE={cfg.get('qtype', 'bi')}")
print(f"CORRECTION={cfg.get('test_correction', 'Bon')}")
print(f"INDEX={cfg.get('index', 'early')}")
print(f"OUTGROUP={cfg.get('outgroup', 'OUT')}")
EOF
)"

echo "alpha=$ALPHA beta=$BETA qtype=$QTYPE correction=$CORRECTION index=$INDEX outgroup=$OUTGROUP"

Rscript ectoblob/ectoblob.R \
    --genes "$SCRATCH/genes.nwk" \
    --tree "$SCRATCH/start.nwk" \
    --out "$SCRATCH/tob.tre" \
    --alpha "$ALPHA" --beta "$BETA" \
    --qtype "$QTYPE" --correction "$CORRECTION" \
    --index "$INDEX" --outgroup "$OUTGROUP"

echo "=== output ==="
head -c 200 "$SCRATCH/tob.tre"
echo
echo "=== test PASSED: output kept in $SCRATCH ==="
