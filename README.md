# Local Sequence Context Influence on Rates of Substitution in the 1000G Deep Sequencing Data

In this repository one can find the scripts we used to analyize the distribution of nucleotides flanking single nucleotide variants (SNVs) in the 1000 Genomes Project Deep Sequencing Data. In particular, we characterize the distributions of nucleotides flanking singletons to make inferences on processes active in modern human populations.

We assume that the vcf files have already been filtered to only include the 2,504 unrelated individuals in the 1000G dataset; the list of ids for these samples is in the file `reference_data/sample_ids.txt`

## Getting GC content in windows

```{bash, eval=FALSE}
curl -s http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes | \
  grep -Ev "_|X|Y|M" > reference_data/hg38.chrom.sizes
  
bedtools makewindows -g reference_data/hg38.chrom.sizes -w 1000000 | grep -Ev "_|X|Y|M" | sort -k 1,1V -k2,2n > reference_data/genome.1000kb.sorted.bed

bedtools makewindows -g reference_data/hg38.chrom.sizes -w 5000000 | grep -Ev "_|X|Y|M" | sort -k 1,1V -k2,2n > reference_data/genome.5000kb.sorted.bed

bedtools makewindows -g reference_data/hg38.chrom.sizes -w 100000 | grep -Ev "_|X|Y|M" | sort -k 1,1V -k2,2n > reference_data/genome.100kb.sorted.bed

bedtools makewindows -g reference_data/hg38.chrom.sizes -w 10000 | grep -Ev "_|X|Y|M" | sort -k 1,1V -k2,2n > reference_data/genome.10kb.sorted.bed

bedtools makewindows -g reference_data/hg38.chrom.sizes -w 1000 | grep -Ev "_|X|Y|M" | sort -k 1,1V -k2,2n > reference_data/genome.1kb.sorted.bed

# GC content
bedtools nuc -fi reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -bed reference_data/genome.1kb.sorted.bed > reference_data/gc1kb.bed

bedtools nuc -fi reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -bed reference_data/genome.10kb.sorted.bed > reference_data/gc10kb.bed

bedtools nuc -fi reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -bed reference_data/genome.100kb.sorted.bed > reference_data/gc100kb.bed
```

# Step 1: Generate singleton files

Here we'll generate singleton files for each of the five super-populations, along with a singleton file for all 2,504 unrelated samples in the 1000G sample. We do this using a chain of calls to `bcftools view` and `vcftools --singletons`. Batch scripts to perform these operations are in the `src` directory in the scripts named `step1_singletons_{POP}.sh`.

# Step 2: Annotate Singletons

For each singleton observation, we pull the 21-mer motif centered at the site of the singleton. For this, we'll need a copy of the human reference genome hg38:

```{bash}
wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gzip -d > reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

The script which appends each singleton with its motif is `src/append_motif.py`, and the batch scripts `step2_annotate_{POP}.sh` will submit slurm jobs for each chromosome for the given population. The output files will be headerless csvs with the following columns:

1. Chromosome
2. Position
3. Original Motif
4. Simple subtype (REF>ALT)
5. ALT
6. Sample ID
7. REF
8. Full Motif (motif and its reverse complement)
9. Condensed sub-type

# Step 3: Sample Control Observations

The scripts to sample 5 controls per singleton are `src/control_sample_at.py` and `src/control_sample_gc.py`. Scripts to submit batch jobs to slurm are available in files with the prefix `step3_sample` (also in the src directory).

NOTE: if using a reference genome other than hg38, you might need to change the variable ref_prefix in the main function definition of the two scripts to match the chromosome names in the fasta file (i.e. in hg37, chromosomes are named only by their number, where as in hg38 each chromosome is prefixed with "chr", e.g. chr1, chr2, ...)

# Step 4: Generate Per-Subtype Files

In the above steps, we've generated files per chromosome. In this step we'll compile singletons and controls across all chromosomes.

Within each subdirectory of `output/singletons`, run the following commands (note: this will take a few minutes to run):

```{bash}
# Generate per-subtype singleton files

awk -F, '{if($9 == "AT_CG")print(substr($8,1,21))}' chr*_annotated.csv > AT_CG.txt
awk -F, '{if($9 == "AT_GC")print(substr($8,1,21))}' chr*_annotated.csv > AT_GC.txt
awk -F, '{if($9 == "AT_TA")print(substr($8,1,21))}' chr*_annotated.csv > AT_TA.txt
awk -F, '{if($9 == "GC_AT")print(substr($8,1,21))}' chr*_annotated.csv > GC_AT.txt
awk -F, '{if($9 == "GC_TA")print(substr($8,1,21))}' chr*_annotated.csv > GC_TA.txt
awk -F, '{if($9 == "GC_CG")print(substr($8,1,21))}' chr*_annotated.csv > GC_CG.txt
awk -F, '{if($9 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_annotated.csv > cpg_GC_AT.txt
awk -F, '{if($9 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_annotated.csv > cpg_GC_TA.txt
awk -F, '{if($9 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_annotated.csv> cpg_GC_CG.txt
```

Within each subdirectory of  `output/controls`, run the following commands to yield per-subtype files:

```{bash}
awk -F, '{if($4 == "AT_CG")print(substr($8,1,21))}' chr*_at.csv > AT_CG.txt
awk -F, '{if($4 == "AT_GC")print(substr($8,1,21))}' chr*_at.csv > AT_GC.txt
awk -F, '{if($4 == "AT_TA")print(substr($8,1,21))}' chr*_at.csv > AT_TA.txt
awk -F, '{if($4 == "GC_AT")print(substr($8,1,21))}' chr*_gc.csv > GC_AT.txt
awk -F, '{if($4 == "GC_TA")print(substr($8,1,21))}' chr*_gc.csv > GC_TA.txt
awk -F, '{if($4 == "GC_CG")print(substr($8,1,21))}' chr*_gc.csv > GC_CG.txt
awk -F, '{if($4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc.csv > cpg_GC_AT.txt
awk -F, '{if($4 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_gc.csv > cpg_GC_TA.txt
awk -F, '{if($4 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_gc.csv> cpg_GC_CG.txt
```

## UPDATE (31-Jan-2022)

Run the following in the `output/controls/ALL` subdirectory to generate the control file for the GC_AT subtype when we ignore CpG status when sampling controls:

```{bash}
awk -F, '{if($4 == "GC_AT" || $4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_all_gc.csv > all_GC_AT.txt
awk -F, '{if($4 == "GC_TA" || $4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_all_gc.csv > all_GC_TA.txt
awk -F, '{if($4 == "GC_CG" || $4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_all_gc.csv > all_GC_CG.txt
```

## UPDATE (10-Nov-2022)

Get subtype files for nearest and furthest control

```{bash}
awk -F, '{if($4 == "AT_CG")print(substr($8,1,21))}' chr*_at.csv.min > AT_CG_min.txt
awk -F, '{if($4 == "AT_GC")print(substr($8,1,21))}' chr*_at.csv.min > AT_GC_min.txt
awk -F, '{if($4 == "AT_TA")print(substr($8,1,21))}' chr*_at.csv.min > AT_TA_min.txt
awk -F, '{if($4 == "GC_AT")print(substr($8,1,21))}' chr*_gc.csv.min > GC_AT_min.txt
awk -F, '{if($4 == "GC_TA")print(substr($8,1,21))}' chr*_gc.csv.min > GC_TA_min.txt
awk -F, '{if($4 == "GC_CG")print(substr($8,1,21))}' chr*_gc.csv.min > GC_CG_min.txt
awk -F, '{if($4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc.csv.min > cpg_GC_AT_min.txt
awk -F, '{if($4 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_gc.csv.min > cpg_GC_TA_min.txt
awk -F, '{if($4 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_gc.csv.min > cpg_GC_CG_min.txt

awk -F, '{if($4 == "AT_CG")print(substr($8,1,21))}' chr*_at.csv.max > AT_CG_max.txt
awk -F, '{if($4 == "AT_GC")print(substr($8,1,21))}' chr*_at.csv.max > AT_GC_max.txt
awk -F, '{if($4 == "AT_TA")print(substr($8,1,21))}' chr*_at.csv.max > AT_TA_max.txt
awk -F, '{if($4 == "GC_AT")print(substr($8,1,21))}' chr*_gc.csv.max > GC_AT_max.txt
awk -F, '{if($4 == "GC_TA")print(substr($8,1,21))}' chr*_gc.csv.max > GC_TA_max.txt
awk -F, '{if($4 == "GC_CG")print(substr($8,1,21))}' chr*_gc.csv.max > GC_CG_max.txt
awk -F, '{if($4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc.csv.max > cpg_GC_AT_max.txt
awk -F, '{if($4 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_gc.csv.max > cpg_GC_TA_max.txt
awk -F, '{if($4 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_gc.csv.max > cpg_GC_CG_max.txt
```

## UPDATE (15 Nov 2022)

Get position files

### Controls

```{bash}
for i in `seq 1 22`; do
echo "AT_CG..."
awk -F, '{if($4 == "AT_CG")print($9)}' "chr${i}_at.csv" >> pos_files/"AT_CG_${i}.txt"
echo "AT_GC..."
awk -F, '{if($4 == "AT_GC")print($9)}' "chr${i}_at.csv" >> pos_files/"AT_GC_${i}.txt"
echo "AT_TA..."
awk -F, '{if($4 == "AT_TA")print($9)}' "chr${i}_at.csv" >> pos_files/"AT_TA_${i}.txt"
done

for i in `seq 1 22`; do
echo "GC_AT..."
awk -F, '{if($4 == "GC_AT")print($9)}' "chr${i}_gc.csv" >> pos_files/"GC_AT_${i}.txt"
echo "GC_TA..."
awk -F, '{if($4 == "GC_TA")print($9)}' "chr${i}_gc.csv" >> pos_files/"GC_TA_${i}.txt"
echo "GC_CG..."
awk -F, '{if($4 == "GC_CG")print($9)}' "chr${i}_gc.csv" >> pos_files/"GC_CG_${i}.txt"
done

for i in `seq 1 22`; do
echo "GC_AT..."
awk -F, '{if($4 == "cpg_GC_AT")print($9)}' "chr${i}_gc.csv" >> pos_files/"cpg_GC_AT_${i}.txt"
echo "GC_TA..."
awk -F, '{if($4 == "cpg_GC_TA")print($9)}' "chr${i}_gc.csv" >> pos_files/"cpg_GC_TA_${i}.txt"
echo "GC_CG..."
awk -F, '{if($4 == "cpg_GC_CG")print($9)}' "chr${i}_gc.csv" >> pos_files/"cpg_GC_CG_${i}.txt"
done
```

### Singletons

```{bash}
for i in `seq 1 22`; do
echo "AT_CG..."
awk -F, '{if($9 == "AT_CG")print($2)}' "chr${i}_annotated.csv" >> pos_files/"AT_CG_${i}.txt"
echo "AT_GC..."
awk -F, '{if($9 == "AT_GC")print($2)}' "chr${i}_annotated.csv" >> pos_files/"AT_GC_${i}.txt"
echo "AT_TA..."
awk -F, '{if($9 == "AT_TA")print($2)}' "chr${i}_annotated.csv" >> pos_files/"AT_TA_${i}.txt"
done

for i in `seq 1 22`; do
echo "GC_AT..."
awk -F, '{if($9 == "GC_AT")print($2)}' "chr${i}_annotated.csv" >> pos_files/"GC_AT_${i}.txt"
echo "GC_TA..."
awk -F, '{if($9 == "GC_TA")print($2)}' "chr${i}_annotated.csv" >> pos_files/"GC_TA_${i}.txt"
echo "GC_CG..."
awk -F, '{if($9 == "GC_CG")print($2)}' "chr${i}_annotated.csv" >> pos_files/"GC_CG_${i}.txt"
done

for i in `seq 1 22`; do
echo "GC_AT..."
awk -F, '{if($9 == "cpg_GC_AT")print($2)}' "chr${i}_annotated.csv" >> pos_files/"cpg_GC_AT_${i}.txt"
echo "GC_TA..."
awk -F, '{if($9 == "cpg_GC_TA")print($2)}' "chr${i}_annotated.csv" >> pos_files/"cpg_GC_TA_${i}.txt"
echo "GC_CG..."
awk -F, '{if($9 == "cpg_GC_CG")print($2)}' "chr${i}_annotated.csv" >> pos_files/"cpg_GC_CG_${i}.txt"
done
```

# Step 5: Genome-wide Background Rates - Single Position Models

The scripts to generate counts based on the reference genome are `gw_1_count3cats.py` and `gw_1_count_6cats.py`. Batch scripts to submit jobs to slurm are `step5_gw_1_count_3cat_batch.sh` and `step5_gw_1_count_6cat_batch.sh`.

After this, run the script `step5_combine_chromosomes.R` to combine the per-chromosome per-reference files into per-reference files.

# Step 6: Genome-wide Background Rates - Two Postion Models

The script to generate counts based on the reference genome is `gw_2_count_6cats.py`, with batch script to submit jobs to slurm `step6_gw_2_6cats.sh`.

After this is run, run the R script `step6_combine_chrom.R`. Then, within the `output/gw_2_count/6st` directory run:

```{bash}
for file in non_*; do mv "$file" "${file#non_}";done;
```

Finally, run the script `step6_3cats.R` to condense sub-categories.


# Addendums

## 24 Jan 2022

* Added scripts to sample GC_NN controls that do not account for CpG/non-CpG status (`src/control_sample_all_gc.py`)
