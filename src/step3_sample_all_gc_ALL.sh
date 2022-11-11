#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --time=01:30:00
#SBATCH --job-name=ALL
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/sas-%a.err
#SBATCH -o slurm/sas-%a.out

POP='ALL'

echo "GC ${SLURM_ARRAY_TASK_ID} $POP"
python /net/snowwhite/home/beckandy/research/1000G_LSCI/src/control_sample_all_gc.py \
-s /net/snowwhite/home/beckandy/research/1000G_LSCI/output/singletons/$POP/chr${SLURM_ARRAY_TASK_ID}_annotated.csv \
-f /net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
-o /net/snowwhite/home/beckandy/research/1000G_LSCI/output/controls/$POP/chr${SLURM_ARRAY_TASK_ID}_gc_all.csv \
-n 5 ${SLURM_ARRAY_TASK_ID}
