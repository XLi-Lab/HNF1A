# this r script ...


# load libraries
library(readr)
library(dplyr)

# read the MA0046.1 FIMO output
GSM6248576_FIMO_MA0046.1 <- read_tsv(
  "~/Downloads/GSE206240_RAW/HNF1A/HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_MA0046.1.tsv",
  comment = "#", show_col_types = FALSE)

# number of peaks containing at least one HNF1A motif
n_motif_peaks <- GSM6248576_FIMO_MA0046.1 %>%
  distinct(sequence_name) %>%
  nrow()

n_motif_peaks

# as a fraction of all peaks
n_total_peaks <- nrow(GSM6248576)
cat("Peaks with motif (MA0046.1):", n_motif_peaks, "/", n_total_peaks,
    sprintf("(%.1f%%)\n", 100 * n_motif_peaks / n_total_peaks))