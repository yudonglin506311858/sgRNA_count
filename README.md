# sgRNA_count
实现mageck的定量，不过是用R语言。

#SOFTWARE
# Create a new conda environment
conda create -n crispr_analysis r-base r-essentials r-biocmanager bioconductor-shortread bioconductor-biostrings r-dplyr r-stringr r-optparse r-purrr

# Activate the environment
conda activate crispr_analysis


# example
Rscript sgRNA_count.R \
  --fastq ctr-1_R1.fq.gz,ctr-2_R1.fq.gz,ctr-3_R1.fq.gz \
  --names ctr-1,ctr-2,ctr-3 \
  --ncpu 40 \
  --library library.csv \
  --output sgRNA_counts
![image](https://github.com/user-attachments/assets/5e57b362-3354-4e24-9507-77eb09f4975f)
![image](https://github.com/user-attachments/assets/5e57b362-3354-4e24-9507-77eb09f4975f)
