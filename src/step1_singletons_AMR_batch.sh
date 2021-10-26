#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=700MB
#SBATCH --ntasks=1
#SBATCH --time 04:00:00
#SBATCH --job-name=step1AMR
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step1_amr-%J.err
#SBATCH -o slurm/step1_amr-%J.out

bcftools view --types 'snps' -f 'FILTER="PASS"' -O u /net/topmed8/working/call_sets/freeze9/release/subset/1000g/minDP10.minAC1.PASS.phased/freeze9.1000g.chr${SLURM_ARRAY_TASK_ID}.filtered.gtonly.minDP10.minAC1.PASS.phased.vcf.gz |\
bcftools view -S /net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/AMR_IDs.txt -x -O u |\
bcftools view -c 1 -C 1 -O v |\
vcftools --vcf - --singletons -c > /net/snowwhite/home/beckandy/research/1000G_LSCI/output/singletons/AMR/chr${SLURM_ARRAY_TASK_ID}.txt 
