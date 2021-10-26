#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=8000MB
#SBATCH --ntasks=1
#SBATCH --time 06:00:00
#SBATCH --job-name=s_AFR
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step2_afr-%J.err
#SBATCH -o slurm/step2_afr-%J.out

python append_motif.py -c ${SLURM_ARRAY_TASK_ID} \
-s /net/snowwhite/home/beckandy/research/1000G_LSCI/output/singletons/AFR/chr${SLURM_ARRAY_TASK_ID}.txt \
-o /net/snowwhite/home/beckandy/research/1000G_LSCI/output/singletons/AFR/chr${SLURM_ARRAY_TASK_ID}_annotated.csv
