#!/bin/bash
#SBATCH --job-name=IIWeek
#SBATCH --array=0-5
#SBATCH --ntasks-per-node=24
#SBATCH --mem=31744
#SBATCH --nodes=1
#SBATCH --time=60:00:00
python Binder.py ./../data/input/pairs_25012019122047.txt $SLURM_ARRAY_TASK_ID 2>&1 | tee -a log_25012019122047.txt