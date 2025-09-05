#!/bin/bash
#SBATCH --job-name=n300xn300finalproj
#SBATCH --partition day,scavenge
#SBATCH --requeue
#SBATCH --mem=2G
#SBATCH --mail-type=none
#SBATCH --array=1-600
#SBATCH -t 6:00:00

module load R
Rscript ./model_sim.R $SLURM_ARRAY_TASK_ID  ${is_diff} ${beta3} ${beta4} ${con} ${beta5} ${alph}