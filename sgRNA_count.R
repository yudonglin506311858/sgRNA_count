#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ShortRead)
  library(parallel)
  library(stringr)
  library(dplyr)
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("-f", "--fastq"), type="character", default=NULL,
              help="Comma-separated list of FASTQ files", metavar="character"),
  make_option(c("-n", "--names"), type="character", default=NULL,
              help="Comma-separated sample names corresponding to FASTQ files", metavar="character"),
  make_option(c("-c", "--ncpu"), type="integer", default=10,
              help="Number of CPU cores to use [default= %default]", metavar="integer"),
  make_option(c("-l", "--library"), type="character", default=NULL,
              help="sgRNA library CSV file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="sgRNA_counts",
              help="Output prefix [default= %default]", metavar="character")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$fastq) | is.null(opt$library)) {
  print_help(opt_parser)
  stop("Both --fastq and --library arguments must be supplied", call.=FALSE)
}

# Record start time
start_time <- Sys.time()

# 1. Load sgRNA library
message("Loading sgRNA library: ", opt$library)
gRNA_library <- read.csv(opt$library, stringsAsFactors = FALSE)
target_sequences <- unique(gRNA_library$gRNA.sequence)

# 2. Process input files and sample names
fastq_files <- unlist(strsplit(opt$fastq, ","))

if (!is.null(opt$names)) {
  sample_names <- unlist(strsplit(opt$names, ","))
  if (length(fastq_files) != length(sample_names)) {
    stop("Number of FASTQ files must match number of sample names", call.=FALSE)
  }
} else {
  # Generate sample names from filenames if not provided
  sample_names <- tools::file_path_sans_ext(basename(fastq_files))
  sample_names <- gsub("_R1\\.fq\\.gz$", "", sample_names)
}

# 3. Initialize final result table
final_result <- gRNA_library %>%
  select(id, gRNA.sequence, Gene)

# 4. Process each FASTQ file
for (i in seq_along(fastq_files)) {
  fastq_file <- fastq_files[i]
  sample_name <- sample_names[i]
  
  message("\nProcessing sample ", i, " of ", length(fastq_files), ": ", sample_name)
  message("  FASTQ file: ", fastq_file)
  
  # Read FASTQ
  message("  Reading FASTQ file...")
  fq <- readFastq(fastq_file)
  reads <- as.character(sread(fq))
  
  # Parallel counting
  message("  Counting sgRNA matches using ", opt$ncpu, " CPU cores...")
  cl <- makeCluster(opt$ncpu, type = "PSOCK")
  clusterExport(cl, varlist = c("reads", "str_detect", "fixed"), envir = environment())
  
  count_list <- parLapply(cl, target_sequences, function(seq) {
    sum(str_detect(reads, fixed(seq)))
  })
  
  stopCluster(cl)
  
  # Create result data frame
  result <- data.frame(
    gRNA.sequence = target_sequences,
    count = unlist(count_list)
  )
  
  # Add to final result
  final_result <- final_result %>%
    left_join(result, by = "gRNA.sequence") %>%
    rename(!!sample_name := count)
  
  # Save individual result
  individual_result <- gRNA_library %>%
    left_join(result, by = "gRNA.sequence") %>%
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    select(id, gRNA.sequence, Gene, count) %>%
    rename(!!sample_name := count)
  
  indiv_file <- paste0(opt$output, "_", sample_name, ".csv")
  message("  Writing individual results to: ", indiv_file)
  write.csv(individual_result, indiv_file, row.names = FALSE)
}

# 5. Final processing and output
message("\nFinalizing results...")
final_result[is.na(final_result)] <- 0

merged_file <- paste0(opt$output, "_merged_all_samples.csv")
message("Writing merged results to: ", merged_file)
write.csv(final_result, merged_file, row.names = FALSE)

# 6. Print runtime
end_time <- Sys.time()
runtime <- end_time - start_time
message("\nAnalysis completed!")
message("Total runtime: ", round(runtime, 2), " ", units(runtime))
