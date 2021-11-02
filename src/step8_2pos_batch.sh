#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=1500MB
#SBATCH --ntasks=1
#SBATCH --time 12:00:00
#SBATCH --job-name=1000G
#SBATCH --partition=nomosix
#SBATCH --array=1-54
#SBATCH --requeue
#SBATCH -e slurm/step8-%J.err
#SBATCH -o slurm/step8-%J.out


Rscript /net/snowwhite/home/beckandy/research/1000G_LSCI/src/2_pos_model_data.R ${SLURM_ARRAY_TASK_ID}
