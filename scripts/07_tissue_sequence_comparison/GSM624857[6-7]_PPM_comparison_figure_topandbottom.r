# this script is builds tissue-specific HNF1A sequence logos 
# and tests per-position base composition between islets and HepG2

# load libraries
library(readr)
library(Biostrings)
library(ggseqlogo)   # for logos
library(dplyr)
library(gridExtra)

GSM6248576.7_Combined_FIMO <- read_csv("GSM624857[6-7]_Outputs/GSM624857[6-7]_Tissue_Sequence_Comparison.csv", show_col_types = FALSE)

# pull sequences for each tissue-specific group
islet_specific_seqs <- GSM6248576.7_Combined_FIMO %>%
  filter(target_group == "Islet-Specific") %>%
  pull(matched_sequence)

hep_specific_seqs <- GSM6248576.7_Combined_FIMO %>%
  filter(target_group == "HepG2-Specific") %>%
  pull(matched_sequence)

# how many sequences in each?
cat("Islet-specific motif hits:", length(islet_specific_seqs), "\n")
cat("HepG2-specific motif hits:", length(hep_specific_seqs), "\n")

# confirm all sequences are the same length (they should be — FIMO uses fixed motif width)
cat("Islet sequence lengths:", paste(unique(nchar(islet_specific_seqs)), collapse = ", "), "\n")
cat("HepG2 sequence lengths:", paste(unique(nchar(hep_specific_seqs)),  collapse = ", "), "\n")

# build the PPMs
# convert to DNAStringSet (Biostrings format)
islet_dna <- DNAStringSet(islet_specific_seqs)
hep_dna   <- DNAStringSet(hep_specific_seqs)

# build position frequency matrices (PFMs)
islet_pfm <- consensusMatrix(islet_dna, as.prob = FALSE)[c("A","C","G","T"), ]
hep_pfm   <- consensusMatrix(hep_dna,   as.prob = FALSE)[c("A","C","G","T"), ]

# convert to position probability matrices (PPMs), divide by column sums
islet_ppm <- sweep(islet_pfm, 2, colSums(islet_pfm), "/")
hep_ppm   <- sweep(hep_pfm,   2, colSums(hep_pfm),   "/")

# each column should sum to 1
colSums(islet_ppm)
colSums(hep_ppm)

# plot the logos side by side
# ggseqlogo accepts the PPM matrix directly
p_islet <- ggseqlogo(islet_ppm, method = "prob") +
  ggtitle(paste0("Islet-Specific HNF1A Motif (n=", length(islet_specific_seqs), ")")) +
  theme(
    plot.title   = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text    = element_text(size = 12, face = "bold"),
    plot.margin  = margin(t = 5, r = 5, b = 5, l = 5))

p_hep <- ggseqlogo(hep_ppm, method = "prob") +
  ggtitle(paste0("HepG2-Specific HNF1A Motif (n=", length(hep_specific_seqs), ")")) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text  = element_text(size = 12, face = "bold"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

grid.arrange(p_islet, p_hep, nrow = 2)

# statistical testing at each position of the sequence
# per-position base counts (4 x 14), built from the tissue-specific matched sequences (FIMO orientation)
cm_islet <- consensusMatrix(DNAStringSet(islet_specific_seqs), as.prob = FALSE)[c("A","C","G","T"), ]
cm_hep   <- consensusMatrix(DNAStringSet(hep_specific_seqs),   as.prob = FALSE)[c("A","C","G","T"), ]

L <- ncol(cm_islet)

res <- do.call(rbind, lapply(seq_len(L), function(p) {
  tab <- rbind(islet = cm_islet[, p], hepG2 = cm_hep[, p])   # 2 x 4
  tab <- tab[, colSums(tab) > 0, drop = FALSE]               # drop bases absent in both
  ex_ok <- all(suppressWarnings(chisq.test(tab)$expected) >= 5)
  test  <- if (ex_ok) chisq.test(tab) else
    chisq.test(tab, simulate.p.value = TRUE, B = 1e4)
  n <- sum(tab)
  V <- sqrt(unname(test$statistic) / (n * (min(dim(tab)) - 1)))   # Cramér's V
  data.frame(position  = p,
             method    = if (ex_ok) "chisq" else "chisq_MC",
             statistic = round(unname(test$statistic), 2),
             p_value   = signif(test$p.value, 3),
             cramers_v = round(V, 3))}))
res$p_adj <- signif(p.adjust(res$p_value, method = "BH"), 3)
res

per_pos <- function(ppm, tissue) {
  data.frame(
    tissue   = tissue,
    position = seq_len(ncol(ppm)),
    base     = rownames(ppm)[apply(ppm, 2, which.max)],  # dominant base
    prob     = round(apply(ppm, 2, max), 2))}            # its probability

# both tissues side by side, one row per position
islet_tbl <- per_pos(islet_ppm, "Islet")
hep_tbl   <- per_pos(hep_ppm,   "HepG2")
merge(islet_tbl, hep_tbl, by = "position", suffixes = c("_islet", "_hep"))

