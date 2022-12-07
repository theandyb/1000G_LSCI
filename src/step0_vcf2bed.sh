#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=2GB
#SBATCH --ntasks=1
#SBATCH --time 02:00:00
#SBATCH --job-name=step0
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/vcf_loc-%J.err
#SBATCH -o slurm/vcf_loc-%J.out

bcftools view -G --types 'snps' -f 'FILTER="PASS"' -O v /net/topmed8/working/call_sets/freeze9/release/subset/1000g/minDP10.minAC1.PASS.phased/freeze9.1000g.chr${SLURM_ARRAY_TASK_ID}.filtered.gtonly.minDP10.minAC1.PASS.phased.vcf.gz |\
bcftools query -f '%CHROM\t%POS\t%POS\n' | awk '{print($1"\t"$2 - 1"\t"$3)}' >\
/net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/vcfbed/chr${SLURM_ARRAY_TASK_ID}.bed

