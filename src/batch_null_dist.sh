#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --time=01:30:00
#SBATCH --job-name=SAS
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/null_sas_%a.err
#SBATCH -o slurm/null_sas_%a.out

POP='SAS'

echo "$POP ${SLURM_ARRAY_TASK_ID}"
python /net/snowwhite/home/beckandy/research/1000G_LSCI/src/null_dist.py \
-s /net/snowwhite/home/beckandy/research/1000G_LSCI/output/singletons/$POP/chr${SLURM_ARRAY_TASK_ID}_annotated.csv \
-f /net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
-o /net/snowwhite/home/beckandy/research/1000G_LSCI/output/null_dist/$POP/ ${SLURM_ARRAY_TASK_ID}
