# Generate bed file with variant locations

# VCF downloaded from http://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/
# bcftools view -> only SNPs

for i in `seq 1 22`; do
echo $i
bcftools view -v snps homo_sapiens-chr$i.vcf.gz | bcftools query -f '%CHROM\t%POS\n' |\
awk '{print($1"\t"$2-1"\t"$2)}' >> mask.bed
done

# mask-a the fasta
bedtools maskfasta -fi /net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -bed /net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/ensembl/mask.bed \
  -fo /net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/masked.fa -soft
