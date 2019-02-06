#!/bin/bash
#SBATCH --job-name=Control
#SBATCH --array=0-3
#SBATCH --ntasks-per-node=24
#SBATCH --mem=31744
#SBATCH --nodes=1
#SBATCH --time=72:00:00
cd vbind/src
python Binder.py ./../data/input/pairs_04022019220519.txt $SLURM_ARRAY_TASK_ID 2>&1 | tee -a log_04022019220519.txt