#!/bin/bash
#SBATCH --job-name=III_EV_RDR
#SBATCH --array=0-3
#SBATCH --ntasks-per-node=24
#SBATCH --mem=31744
#SBATCH --nodes=1
#SBATCH --time=72:00:00
python Binder.py ./../data/input/pairs_24012019133600.txt $SLURM_ARRAY_TASK_ID 2>&1 | tee -a log_24012019133600.txt