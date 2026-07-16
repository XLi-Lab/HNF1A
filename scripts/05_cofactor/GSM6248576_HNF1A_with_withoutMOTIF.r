# this r script is important before conduct de novo motif analysis.
# creates two FASTA files that will be the input for XSTREME
# one file containing sequences of peaks WITH HNF1A motif, and one file containing sequences of peaks WITHOUT HNF1A motif

# load libraries
library(Biostrings)
library(dplyr)
library(readr)

# paths
GSM6248576_FASTA <- "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_peaks-FASTAs.fa"
GSM6248576_FIMO  <- "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_MA0046.1.tsv"

HNF1A_With    <- "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_HNF1A_withMotif.fa"
HNF1A_Without <- "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_HNF1A_withoutMotif.fa"

# load peak sequences
seqs <- readDNAStringSet(GSM6248576_FASTA)

# load FIMO results
GSM6248576_MA00461_FIMO <- read_tsv(GSM6248576_FIMO, comment = "#", show_col_types = FALSE)

# FIMO column is usually called sequence_name
motif_peak_names <- unique(GSM6248576_MA00461_FIMO$sequence_name)

# split sequences
with_motif <- seqs[names(seqs) %in% motif_peak_names]
without_motif <- seqs[!(names(seqs) %in% motif_peak_names)]

# export
writeXStringSet(with_motif,    HNF1A_With,    format = "fasta")
writeXStringSet(without_motif, HNF1A_Without, format = "fasta")

# quick check
cat("Total peaks:", length(seqs), "\n")
cat("With HNF1A motif:", length(with_motif), "\n")
cat("Without HNF1A motif:", length(without_motif), "\n")

