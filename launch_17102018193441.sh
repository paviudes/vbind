#SBATCH --jobname=IIIweek
#SBATCH --array=0-6
#SBATCH --cpus-per-task=4
#SBATCH --mem=32768
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=120:00:00
python Binder.py pairs_17102018193441.txt $SLURM_ARRAY_TASK_ID 2>&1 | tee -a log.txt