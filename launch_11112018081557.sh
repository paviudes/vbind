#!/bin/bash
#SBATCH --job-name=unfinished
#SBATCH --array=0-4
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=120:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pavithran.sridhar@gmail.com
python Binder.py pairs_11112018081557.txt $SLURM_ARRAY_TASK_ID 2>&1 | tee -a log.txt
