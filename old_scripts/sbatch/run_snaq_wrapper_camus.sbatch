#!/usr/bin/env bash
#SBATCH --time=7-00:00:00
#SBATCH --job-name=run_all_snaq
#SBATCH --partition=tallis
#SBATCH --account=Tallis
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g 
#SBATCH --output=/u/quadros3/research/phylogeny/project_01/experiment_01/1a/testing/logs/run_all_snaq_%j_log.txt

/u/quadros3/research/phylogeny/project_01/experiment_01/1a/scripts/bash/run_snaq_wrapper.sh /u/quadros3/research/phylogeny/project_01/experiment_01/1a/data/camus-dataset/n15/ /u/quadros3/research/phylogeny/project_01/experiment_01/1a/data/inferred-networks/n15/ /u/quadros3/research/phylogeny/project_01/experiment_01/1a/data/inferred-snaq/n15 1 32 iqtree_500.nwk camus_f05_astral_iqtree.nwk