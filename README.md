# Local Sequence Context Influence on Rates of Substitution in the 1000G Deep Sequencing Data

In this repository one can find the scripts we used to analyize the distribution of nucleotides flanking single nucleotide variants (SNVs) in the 1000 Genomes Project Deep Sequencing Data. In particular, we characterize the distributions of nucleotides flanking singletons to make inferences on processes active in modern human populations.

We assume that the vcf files have already been filtered to only include the 2,504 unrelated individuals in the 1000G dataset; the list of ids for these samples is in the file `reference_data/sample_ids.txt`

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