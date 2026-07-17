# DEAD SCRIPT
# this r script is to conduct reverse motif analysis. for example, of the sequences where HNF1A does not have any motifs, are there any motifs for other interaction partners

# load libraries
library(Biostrings)
library(dplyr)
library(readr)

# paths
GSM6248577_FASTA <- "HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_peaks-FASTAs.fa"
fimo_file  <- "HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_MA0046.3.tsv"

HNF1A_With    <- "HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_HNF1A_withMotif.fa"
HNF1A_Without <- "HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_HNF1A_withoutMotif.fa"

# load peak sequences
seqs <- readDNAStringSet(GSM6248577_FASTA)

# load FIMO results
GSM6248577_MA00463_FIMO <- read_tsv(fimo_file, comment = "#", show_col_types = FALSE)

# FIMO column is usually called sequence_name
motif_peak_names <- unique(GSM6248577_MA00463_FIMO$sequence_name)

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