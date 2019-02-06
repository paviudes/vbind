#!/bin/bash
#SBATCH --job-name=IIIWeek
#SBATCH --array=0-3
#SBATCH --ntasks-per-node=24
#SBATCH --mem=31744
#SBATCH --nodes=1
#SBATCH --time=60:00:00
python Binder.py ./../data/input/pairs_26012019121935.txt $SLURM_ARRAY_TASK_ID 2>&1 | tee -a log_26012019121935.txt