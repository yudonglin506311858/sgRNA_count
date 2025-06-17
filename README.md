# sgRNA_count

# example
Rscript sgRNA_count.R \
  --fastq ctr-1_R1.fq.gz,ctr-2_R1.fq.gz,ctr-3_R1.fq.gz \
  --names ctr-1,ctr-2,ctr-3 \
  --ncpu 40 \
  --library library.csv \
  --output sgRNA_counts
