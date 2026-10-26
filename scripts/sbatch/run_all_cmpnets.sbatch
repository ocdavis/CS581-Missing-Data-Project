#!/usr/bin/env bash
#SBATCH --time=7-00:00:00
#SBATCH --job-name=run-calgtprob-all
#SBATCH --partition=cs
#SBATCH --account=cs-warnow
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --output=/u/quadros3/research/phylogeny/project_01/experiment_01/1a/testing/logs/test-calgtprob-all_%j_log.txt

module load java/23
python3  /u/quadros3/research/phylogeny/project_01/experiment_01/1a/scripts/src/run_all_cmpnets.py -cn /u/quadros3/research/phylogeny/project_01/experiment_01/1a/data/inferred-networks/n15 -tn /u/quadros3/research/phylogeny/project_01/experiment_01/1a/data/camus-dataset/n15 -o /u/quadros3/research/phylogeny/project_01/experiment_01/1a/output/n15 -m cluster -p /u/quadros3/research/phylogeny/project_01/experiment_01/1a/PhyloNet.jar