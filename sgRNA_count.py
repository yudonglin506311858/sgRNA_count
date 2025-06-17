#!/usr/bin/env python3
import gzip
import argparse
from multiprocessing import Pool
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Count sgRNA sequences in FASTQ files')
    parser.add_argument('-f', '--fastq', required=True,
                        help='Comma-separated list of FASTQ files')
    parser.add_argument('-l', '--library', required=True,
                        help='sgRNA library CSV file')
    parser.add_argument('-n', '--names', default=None,
                        help='Comma-separated sample names')
    parser.add_argument('-c', '--ncpu', type=int, default=10,
                        help='Number of CPU cores to use')
    parser.add_argument('-o', '--output', default='sgRNA_counts',
                        help='Output prefix')
    return parser.parse_args()

def load_library(library_file):
    return pd.read_csv(library_file)

def count_sequences(args):
    seq, fastq_file = args
    count = 0
    with gzip.open(fastq_file, 'rt') as f:
        for record in SeqIO.parse(f, 'fastq'):
            if seq in str(record.seq):
                count += 1
    return (seq, count)

def process_sample(fastq_file, library, ncpu):
    targets = library['gRNA.sequence'].unique()
    
    with Pool(ncpu) as pool:
        results = list(tqdm(pool.imap(count_sequences, 
                                    [(seq, fastq_file) for seq in targets]),
                      total=len(targets),
                      desc=f"Processing {os.path.basename(fastq_file)}"))
    
    return pd.DataFrame(results, columns=['gRNA.sequence', 'count'])

def main():
    args = parse_args()
    
    # Load library
    library = load_library(args.library)
    
    # Process samples
    fastq_files = args.fastq.split(',')
    sample_names = args.names.split(',') if args.names else [
        os.path.basename(f).replace('_R1.fq.gz', '').replace('.fq.gz', '') 
        for f in fastq_files
    ]
    
    merged_counts = library.copy()
    
    for fastq_file, name in zip(fastq_files, sample_names):
        print(f"\nProcessing {name} ({fastq_file})")
        counts = process_sample(fastq_file, library, args.ncpu)
        merged_counts = merged_counts.merge(
            counts, on='gRNA.sequence', how='left'
        ).rename(columns={'count': name})
        
        # Save individual results
        individual = library.merge(counts, on='gRNA.sequence', how='left')
        individual[name] = individual['count'].fillna(0)
        individual[['id', 'gRNA.sequence', 'Gene', name]].to_csv(
            f"{args.output}_{name}.csv", index=False
        )
    
    # Save merged results
    merged_counts.fillna(0, inplace=True)
    merged_counts.to_csv(f"{args.output}_merged_all_samples.csv", index=False)
    print("\nAnalysis completed!")

if __name__ == '__main__':
    main()
