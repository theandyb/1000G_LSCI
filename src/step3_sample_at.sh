#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --time=01:30:00
#SBATCH --job-name=EAS
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step3_ATeas-%a.err
#SBATCH -o slurm/step3_ATeas-%a.out

POP='SAS'

echo "AT ${SLURM_ARRAY_TASK_ID}"
python /net/snowwhite/home/beckandy/research/1000G_LSCI/src/control_sample_at.py \
-s /net/snowwhite/home/beckandy/research/1000G_LSCI/output/singletons/$POP/chr${SLURM_ARRAY_TASK_ID}_annotated.csv \
-f /net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/grch38_1kgp_var_mask.fa \
-o /net/snowwhite/home/beckandy/research/1000G_LSCI/output/controls/$POP/chr${SLURM_ARRAY_TASK_ID}_at.csv \
-n 5 ${SLURM_ARRAY_TASK_ID}
