# this r script extracts the hg38 genomic sequence of each HNF1A peak and writes them to a FASTA file (input for FIMO)

# load libraries
library(readr)
library(dplyr)
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

# set working directory
setwd("~/Downloads/GSE206240_RAW/HNF1A")

# input CSV data with correct columns
csv_path <- "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_peakcentric.csv"
GSM6248576 <- read_csv(csv_path, show_col_types = FALSE)

# convert to GRanges using full peak starts and ends 
gr <- GRanges(
  seqnames = GSM6248576$peak_chr,
  ranges   = IRanges(start = GSM6248576$peak_start, end = GSM6248576$peak_end))

# extract sequences from hg38
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)

# name each sequence by peak_id for later merge
names(seqs) <- gsub("[:]", "_", GSM6248576$peak_id)

# write out the FASTA file for individual peak binding sequences
writeXStringSet(seqs, filepath = "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_peaks-FASTAs.fa", format = "fasta")
