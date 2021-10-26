# Local Sequence Context Influence on Rates of Substitution in the 1000G Deep Sequencing Data

In this repository one can find the scripts we used to analyize the distribution of nucleotides flanking single nucleotide variants (SNVs) in the 1000 Genomes Project Deep Sequencing Data. In particular, we characterize the distributions of nucleotides flanking singletons to make inferences on processes active in modern human populations.

We assume that the vcf files have already been filtered to only include the 2,504 unrelated individuals in the 1000G dataset; the list of ids for these samples is in the file `reference_data/sample_ids.txt`

# Step 1: Generate singleton files

Here we'll generate singleton files for each of the five super-populations, along with a singleton file for all 2,504 unrelated samples in the 1000G sample. We do this using a chain of calls to `bcftools view` and `vcftools --singletons`. Batch scripts to perform these operations are in the `src` directory in the scripts named `step1_singletons_{POP}.sh`.
