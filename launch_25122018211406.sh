#!/bin/bash
#SBATCH --job-name=test
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1024
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
python Binder.py pairs_25122018211406.txt $SLURM_ARRAY_TASK_ID 2>&1 | tee -a log.txt