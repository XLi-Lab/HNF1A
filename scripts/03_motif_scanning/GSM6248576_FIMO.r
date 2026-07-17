# this r script parses the islet FIMO MA0046.1 output, reconstructs orientation-aware
# flanking sequences, collapses both-strand/overlapping hits into distinct motif sites,
# joins counts onto the full peak table and labels motif presence

# load libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(Biostrings)
library(readr)
library(GenomicRanges)

# define consensus
consensus <- "GGTTAATNATTAAC"

# read peak centric .csv
GSM6248576 <- read_csv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_peakcentric.csv", show_col_types = FALSE)

# read .tsv file from FIMO output
GSM6248576_FIMO <- read_tsv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_MA0046.1.tsv", comment = "#", show_col_types = FALSE)

# rename problematic columns (due to - instead of _)
GSM6248576_FIMO <- GSM6248576_FIMO %>%
  rename(p_value = `p-value`, q_value = `q-value`)

# read the original FASTA file to pull seqeuences from
GSM6248576_FASTA <- readDNAStringSet("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_peaks-FASTAs.fa")

# extract a window of 3 bases before and after each match, reverse-complementing
# minus-strand hits so flank_seq is always in motif orientation
GSM6248576_FIMO <- GSM6248576_FIMO %>%
  rowwise() %>%
  mutate(
    full_seq = as.character(GSM6248576_FASTA[[sequence_name]]),
    raw_seq  = substr(full_seq, pmax(1, start - 3), pmin(nchar(full_seq), stop + 3)),
    flank_seq = ifelse(
      strand == "-",
      as.character(reverseComplement(DNAString(raw_seq))),
      raw_seq)) %>%
  ungroup()

# naive per-peak hit count (kept for comparison)
GSM6248576_FIMO <- GSM6248576_FIMO %>%
  add_count(sequence_name, name = "motif_hits_per_peak")

# distinct-site count
# collapse both-strand / overlapping FIMO hits
distinct_counts <- GSM6248576_FIMO %>%
  group_by(sequence_name) %>%
  group_modify(~{
    gr <- GRanges(seqnames = "x",
                  ranges = IRanges(start = .x$start, end = .x$stop))
    tibble(distinct_motifs = length(reduce(gr)))}) %>%
  ungroup()

# compare naive vs distinct-site counts
table(GSM6248576_FIMO %>% distinct(sequence_name, motif_hits_per_peak) %>% pull(motif_hits_per_peak))
table(distinct_counts$distinct_motifs)

# join distinct-site counts onto the full peak list (peaks with no hit -> 0)
all_peaks_with_counts <- GSM6248576 %>%
  mutate(sequence_name = sub(":", "_", peak_id)) %>%
  left_join(distinct_counts, by = "sequence_name") %>%
  mutate(distinct_motifs = replace_na(distinct_motifs, 0))

# work out motif-count distribution
motif_distribution <- all_peaks_with_counts %>%
  count(distinct_motifs)

# plot motif distribution bar chart
ggplot(motif_distribution, aes(x = distinct_motifs, y = n)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  xlab("Number of Distinct Motif Sites Per Peak") +
  ylab("Number of Peaks")

# adding motif grouping (for subsequent plotting)
all_peaks_with_counts <- all_peaks_with_counts %>%
  mutate(
    motif_group = ifelse(distinct_motifs > 0, "Motif", "No Motif"),
    motif_group = factor(motif_group, levels = c("No Motif", "Motif")))

# save output
write.csv(all_peaks_with_counts, "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_FIMO_AllPeaks.csv", row.names = FALSE)

