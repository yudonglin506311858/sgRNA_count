# sgRNA_count
用R语言和Python实现mageck的sgRNA定量。

A high-performance R script that rapidly quantifies CRISPR sgRNA abundances from sequencing data by parallel pattern matching against a reference library, with command-line interface for batch processing."

Key features highlighted:

CRISPR-specific: Designed for sgRNA counting

Parallelized: Uses multicore processing (40+ CPUs)

CLI-enabled: Simple command-line execution

Batch-compatible: Processes multiple samples simultaneously

# Create a new conda environment
conda create -n crispr_analysis r-base r-essentials r-biocmanager bioconductor-shortread bioconductor-biostrings r-dplyr r-stringr r-optparse r-purrr

conda install biopython pandas tqdm
conda  install pyfastx

# Activate the environment
conda activate crispr_analysis

# Parameter	Shorthand	Required	Default	Description

--fastq	-f	Yes	-	Comma-separated list of gzipped FASTQ files (e.g., sample1.fq.gz,sample2.fq.gz)

--library	-l	Yes	-	CSV file containing sgRNA library (columns: id, gRNA.sequence, Gene)

--names	-n	No	Auto-generated	Comma-separated sample names matching FASTQ files (e.g., ctrl1,ctrl2,exp1). If omitted, names will be derived from filenames.

--ncpu	-c	No	10	Number of CPU cores for parallel processing (recommend 20-50 for large datasets)

--output	-o	No	sgRNA_counts	Output filename prefix (e.g., specifying --output results creates results_merged_all_samples.csv)

--seed	-s	No	123	Random seed for reproducibility (used in subsampling if implemented)


![image](https://github.com/user-attachments/assets/04359169-dc18-4059-8686-3886fbce98be)


# example
Rscript sgRNA_count.R \
  --fastq ctr-1_R1.fq.gz,ctr-2_R1.fq.gz,ctr-3_R1.fq.gz \
  --names ctr-1,ctr-2,ctr-3 \
  --ncpu 40 \
  --library library.csv \
  --output sgRNA_counts


python sgRNA_count.py  \
  --fastq ctr-1_R1.fq.gz,ctr-2_R1.fq.gz,ctr-3_R1.fq.gz \
  --names ctr-1,ctr-2,ctr-3 \
  --ncpu 40 \
  --library library.csv \
  --output sgRNA_counts
![image](https://github.com/user-attachments/assets/5e57b362-3354-4e24-9507-77eb09f4975f)
