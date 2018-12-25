#!/bin/bash
#SBATCH --job-name=IIIweekCtrl
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8192
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
module spider scipy-stack/2018b;python Binder.py pairs_24122018230005.txt $SLURM_ARRAY_TASK_ID 2>&1 | tee -a log.txt
