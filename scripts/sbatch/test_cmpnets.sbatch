#!/usr/bin/env bash
#SBATCH --time=7-00:00:00
#SBATCH --job-name=test-cmpnets
#SBATCH --partition=cs
#SBATCH --account=cs-warnow
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=quadros3@illinois.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --output=/u/quadros3/research/phylogeny/project_01/experiment_01/1a/testing/logs/test-cmpnets_%j_log.txt

module load java/23
python3 /u/quadros3/research/phylogeny/project_01/experiment_01/1a/scripts/src/run_cmpnets.py -tn /u/quadros3/research/phylogeny/project_01/experiment_01/1a/testing/true_net.nwk -cn /u/quadros3/research/phylogeny/project_01/experiment_01/1a/testing/camus_f05_astral_iqtree.nwk -m cluster -p /u/quadros3/research/phylogeny/project_01/experiment_01/1a/PhyloNet.jar -l /u/quadros3/research/phylogeny/project_01/experiment_01/1a/testing/testing1.log 