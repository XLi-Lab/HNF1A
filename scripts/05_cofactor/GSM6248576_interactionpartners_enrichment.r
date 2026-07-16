# this r script is about seeing if other HNF1A interaction partners bind to the peaks with no HNF1A motif,
# is there enrichment of partner factors in these peaks that could account for the presence of HNF1A?
# script compares each candidate partner-TF motif frequency between HNF1A motif-containing and motif-lacking islet peaks

# load libraries
library(dplyr)
library(tidyr)
library(readr)
library(purrr)

# read the per-peak table and define the with/without HNF1A-motif groups
all_peaks <- read_csv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_FIMO_AllPeaks.csv", show_col_types = FALSE)

# label ALL peaks
peak_groups <- all_peaks %>%
  select(sequence_name, motif_group) %>%
  distinct() %>%
  mutate(group = ifelse(motif_group == "Motif", "HNF1A", "No_HNF1A")) %>%
  select(sequence_name, group)

# load all other TF FIMO files
files <- c(
  "HNF1B_fimo.tsv",
  "HNF4G_fimo.tsv",
  "HNF4A_fimo.tsv",
  "ONECUT1_fimo.tsv",
  "FOXA2_fimo.tsv",
  "FOXA3_fimo.tsv",
  "FOSL1_fimo.tsv",
  "JDP2_fimo.tsv")

fimo_dir <- "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/Partner_TF_FIMO"

fimo_all <- map_dfr(files, function(f) {
  read_tsv(file.path(fimo_dir, f), comment = "#", show_col_types = FALSE) %>%
    mutate(TF = gsub("_fimo.tsv", "", f))})

# one hit per peak per TF
peak_hits <- fimo_all %>%
  distinct(sequence_name, TF)

# attach the HNF1A / No_HNF1A group to each partner-TF hit
peak_hits_grouped <- peak_hits %>%
  left_join(peak_groups, by = "sequence_name")

# safety check, partner-TF peak names must match peak names
if (sum(is.na(peak_hits_grouped$group)) > 0) {
  warning(sum(is.na(peak_hits_grouped$group)),
          " partner-TF hits did not match a peak name and will be dropped — check naming consistency")}

# count peaks carrying each partner-TF motif, per group
summary <- peak_hits_grouped %>%
  filter(!is.na(group)) %>%
  group_by(group, TF) %>%
  summarise(peaks_with_motif = n(), .groups = "drop")

# get total peaks per group
group_sizes <- peak_groups %>%
  group_by(group) %>%
  summarise(total_peaks = n())

summary <- summary %>%
  left_join(group_sizes, by = "group") %>%
  mutate(percent = 100 * peaks_with_motif / total_peaks)

# reshape wide to build a contingency table per TF
summary_wide <- summary %>%
  select(group, TF, peaks_with_motif, total_peaks) %>%
  pivot_wider(
    names_from = group,
    values_from = c(peaks_with_motif, total_peaks))

# Fisher's exact test per TF (motif present/absent x HNF1A / No_HNF1A)
results <- summary_wide %>%
  rowwise() %>%
  mutate(
    # build contingency table
    test = list(fisher.test(matrix(
      c(peaks_with_motif_HNF1A,
        total_peaks_HNF1A - peaks_with_motif_HNF1A,
        peaks_with_motif_No_HNF1A,
        total_peaks_No_HNF1A - peaks_with_motif_No_HNF1A),
      nrow = 2,
      byrow = TRUE))),
    
    p_value = test$p.value,
    odds_ratio = test$estimate) %>%
  select(TF, p_value, odds_ratio)

# BH-adjust across the tested TFs
results$p_adj <- p.adjust(results$p_value, method = "BH")

# note: HNF1B enrichment reflects HNF1A/HNF1B motif-sequence similarity (co-detection),
# not necessarily indirect binding, please interpret with caution (see Discussion).

results <- results %>% arrange(p_adj)
results

