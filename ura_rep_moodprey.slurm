#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=28  # Number of cores - 28 is how many it should be when not on DEV
#SBATCH -t 72:00:00
#SBATCH -p standard
#SBATCH -A politics_dept_cluster

# Run program
module purge
module load gcc R/3.3.0

Rscript ura_replication_moodprey.R ${SLURM_CPUS_PER_TASK}