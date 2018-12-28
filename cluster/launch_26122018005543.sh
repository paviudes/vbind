#!/bin/bash
#SBATCH --job-name=IIIweekCtrl
#SBATCH --array=0-1
#SBATCH --mem=31G
#SBATCH --ntasks-per-node=6
#SBATCH --time=48:00:00
python Binder.py pairs_26122018005543.txt $SLURM_ARRAY_TASK_ID 2>&1 | tee -a log_26122018005543.txt
