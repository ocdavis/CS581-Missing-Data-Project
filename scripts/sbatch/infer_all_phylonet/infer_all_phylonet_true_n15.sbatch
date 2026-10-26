#!/usr/bin/env bash
#SBATCH --time=4:00:00
#SBATCH --job-name=infer_all_phylonet_true_n15
#SBATCH --partition=tallis
#SBATCH --account=Tallis
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g 
#SBATCH --output=/u/quadros3/research/phylogeny/project_01/experiment_01/1a/testing/logs/infer_all_phylonet_true_n15_%j_log.txt

module load java/23
python3  /u/quadros3/research/phylogeny/project_01/experiment_01/1a/scripts/src/infer_all_phylonet.py -c /u/quadros3/research/phylogeny/project_01/experiment_01/1a/scripts/config/infer_all_phylonet/infer_all_phylonet_true_n15.toml
