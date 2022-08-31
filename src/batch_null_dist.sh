#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --time=01:30:00
#SBATCH --job-name=AFR
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/null_afr_%a.err
#SBATCH -o slurm/null_afr_%a.out

POP='AFR'

echo "$POP ${SLURM_ARRAY_TASK_ID}"
python /net/snowwhite/home/beckandy/research/1000G_LSCI/src/null_dist.py \
-s /net/snowwhite/home/beckandy/research/1000G_LSCI/output/singletons/$POP/chr${SLURM_ARRAY_TASK_ID}_annotated.csv \
-f /net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/masked.fa \
-o /net/snowwhite/home/beckandy/research/1000G_LSCI/output/null_dist/$POP/ ${SLURM_ARRAY_TASK_ID}
