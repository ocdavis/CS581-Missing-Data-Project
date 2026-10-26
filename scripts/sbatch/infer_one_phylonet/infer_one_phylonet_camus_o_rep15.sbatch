#!/usr/bin/env bash
#SBATCH --time=7-00:00:00
#SBATCH --job-name=infer_one_phylonet_camus_o_rep15
#SBATCH --partition=cs
#SBATCH --account=cs-warnow
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g 
#SBATCH --output=/u/quadros3/research/phylogeny/project_01/experiment_01/1a/testing/logs/infer_one_phylonet_camus_o_rep15_%j_log.txt

module load java/23
python3  /u/quadros3/research/phylogeny/project_01/experiment_01/1a/scripts/src/infer_one_phylonet.py -c /u/quadros3/research/phylogeny/project_01/experiment_01/1a/scripts/config/infer_one_phylonet/infer_one_phylonet_camus_o_rep15.toml
