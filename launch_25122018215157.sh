#!/bin/bash
#SBATCH --job-name=IIweek
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=8192
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
python Binder.py pairs_25122018215157.txt $SLURM_ARRAY_TASK_ID 2>&1 | tee -a log.txt