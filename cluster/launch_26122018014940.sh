#!/bin/bash
#SBATCH --job-name=IIweekEV
#SBATCH --array=0-1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=31744
#SBATCH --nodes=1
#SBATCH --time=48:00:00
python Binder.py pairs_26122018014940.txt $SLURM_ARRAY_TASK_ID 2>&1 | tee -a log_26122018014940.txt