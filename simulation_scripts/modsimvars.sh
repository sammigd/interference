#!/bin/bash
#SBATCH --job-name=n300xn300finalproj
#SBATCH --partition day,scavenge
#SBATCH --requeue
#SBATCH --mem=2G
#SBATCH --mail-type=none
#SBATCH --array=1-600
#SBATCH -t 6:00:00
#SBATCH --output=/nfs/roberts/project/pi_lf474/sgd37/interference/logs/%A_%a.out
#SBATCH --error=/nfs/roberts/project/pi_lf474/sgd37/interference/logs/%A_%a.err
echo "ngamma is: ${ngamma}"
module load R
Rscript ./model_sim.R $SLURM_ARRAY_TASK_ID  ${is_diff} ${beta3} ${beta4} ${con} ${beta5} ${alph} ${ngamma}