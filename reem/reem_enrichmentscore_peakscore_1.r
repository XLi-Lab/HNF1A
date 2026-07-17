# ============================================================
# Reem × HepG2 sequence matching
# Match Reem's 63 B1H sequences against HepG2 ChIP-seq peak FASTAs,
# pull peak scores, and produce a combined table for plotting.
# ============================================================

library(Biostrings)
library(readr)
library(dplyr)

# ── PATHS ──────────────────────────────────────────────────────────────────────
hepg2_fa_path  <- "/Users/Manu/Downloads/GSE206240_RAW/HNF1A/HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_peaks-FASTAs.fa"
hepg2_csv_path <- "/Users/Manu/Downloads/GSE206240_RAW/HNF1A/HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_FIMO_AllPeaks.csv"

# Islet paths (for reference — islet results already known, but include if you want to rerun)
islet_fa_path  <- "/Users/Manu/Downloads/GSE206240_RAW/HNF1A/HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_peaks-FASTAs.fa"
islet_csv_path <- "/Users/Manu/Downloads/GSE206240_RAW/HNF1A/HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_FIMO_AllPeaks.csv"

# ── REEM'S SEQUENCES ───────────────────────────────────────────────────────────
# Manually entered from the .numbers file (delta_DBD_minus_Lib = enrichment score)
reem_df <- tibble::tribble(
  ~reem_id, ~sequence,                 ~delta,  ~tier,
  0,  "AGTTAATTATTAACCAA",   3.629,  "Tier 1 (bulletproof DBD-specific)",
  2,  "CTCGGTTGATCACTCACTTC", -1.330, "Self-activator (enriched both conditions)",
  3,  "TGAAGTTACTCATTAGCAGG", 1.589,  "Tier 2 (reproducible DBD-specific)",
  4,  "AATATTAATTTTTAACTTCT", -0.099, "No DBD-specific enrichment",
  5,  "TGCATTTAACTATTAACAAT", 3.461,  "Tier 1 (bulletproof DBD-specific)",
  6,  "CATGCAAATTATTGACCTCA", 0.180,  "Self-activator (enriched both conditions)",
  7,  "CTTGTTAATAAGTAATTCCC", -0.131, "No DBD-specific enrichment",
  9,  "CACATTAATTAGTAACTTTT", 0.232,  "No DBD-specific enrichment",
  10, "CTGGGCAAATCATTAACCTC", -1.301, "No DBD-specific enrichment",
  11, "AAAGGTAATGATTAATCTGT", -0.546, "No DBD-specific enrichment",
  12, "TTGGGTTAATCATCACTGGG", 0.506,  "No DBD-specific enrichment",
  13, "AAGGGTAAATCTTTAACGAG", -0.521, "No DBD-specific enrichment",
  14, "GATGTTAATCAGTAACCGTC", -0.185, "No DBD-specific enrichment",
  15, "GCTGGTGGATTATTAACAGA", -0.492, "Self-activator (enriched both conditions)",
  16, "TTTGGTTAGTAATTACTAAA", 2.565,  "Tier 1 (bulletproof DBD-specific)",
  17, "GCTGCCAATCATTAACCTCC", -1.051, "No DBD-specific enrichment",
  18, "GCTGTTAATCATTCATTGGG", 0.062,  "No DBD-specific enrichment",
  19, "TTAAATTAAAGATTAACAAA", 3.339,  "Tier 1 (bulletproof DBD-specific)",
  20, "TTTGTTTAATGATCAACAGT", 0.663,  "No DBD-specific enrichment",
  21, "CCAAGTCAATGTTTAACTGG", -0.716, "No DBD-specific enrichment",
  22, "AGAGAAAATGATTAAGCTTT", -1.131, "No DBD-specific enrichment",
  23, "TAACTTAATATTTAACCTAA", 1.701,  "Tier 2 (reproducible DBD-specific)",
  24, "ATGGTTTAATGATTAACAGC", 0.700,  "Self-activator (enriched both conditions)",
  25, "CTTGTTAATTATTAAACCCT", -0.874, "No DBD-specific enrichment",
  26, "ATTAGACAATCATTAACAAG", 1.684,  "Tier 2 (reproducible DBD-specific)",
  27, "TAAGGATAATTGTTAATGTA", 1.547,  "Tier 2 (reproducible DBD-specific)",
  28, "AATGATTAATAAATAATATG", 0.547,  "No DBD-specific enrichment",
  29, "CACAGTCAAACATTAACTTG", 0.739,  "No DBD-specific enrichment",
  30, "AATGTTAATAGTTAATTTGA", 3.608,  "Tier 1 (bulletproof DBD-specific)",
  31, "CCTGTTTAATAATTACCAGA", 2.933,  "Tier 1 (bulletproof DBD-specific)",
  32, "CAAGTTACTGATTAACCCTA", 0.139,  "No DBD-specific enrichment",
  33, "CTGAATTAATAATGAACAAT", 0.635,  "No DBD-specific enrichment",
  34, "CAAGGTTAACGATTAAATGG", 0.790,  "No DBD-specific enrichment",
  35, "TCTAGTTAATAATCTACAAT", 1.731,  "Tier 2 (reproducible DBD-specific)",
  36, "CATGTTAAATATTAATTTGT", -0.389, "No DBD-specific enrichment",
  37, "GACGTGAATTAGTAACCACC", -1.163, "No DBD-specific enrichment",
  38, "GAAGTTAAAAATTAAGCGGC", 2.519,  "Tier 1 (bulletproof DBD-specific)",
  39, "ACTGTTAATAATTAACACTC", 1.144,  "No DBD-specific enrichment",
  40, "TGTGGTTAATGATTAATCGC", -0.954, "No DBD-specific enrichment",
  41, "GCAGCTAATAATAAACCAGT", -0.818, "No DBD-specific enrichment",
  42, "GCTAGTTAATGATTAGTGAA", 0.631,  "No DBD-specific enrichment",
  43, "GTGGTGGATCATTAACCAGC", 0.179,  "No DBD-specific enrichment",
  44, "ACGGCTTAATGATTAACTAT", 0.421,  "No DBD-specific enrichment",
  45, "CATGTTAATAACAAATCACA", -1.010, "No DBD-specific enrichment",
  46, "TCATGTTAAACATTAACAGC", 1.415,  "No DBD-specific enrichment",
  47, "AGAGGTTGATTTTTAACTAC", -0.683, "Self-activator (enriched both conditions)",
  48, "CAGGGTTAATCATTTCCACG", -0.886, "No DBD-specific enrichment",
  49, "AGGGATAAATTATTAGCAGC", 1.493,  "No DBD-specific enrichment",
  50, "AATGTTAACGATCACCCCAA", 0.701,  "No DBD-specific enrichment",
  51, "AATGGTTAATGTTTAAGCGC", 0.107,  "No DBD-specific enrichment",
  52, "CCCTCTAACCATTAACCACC", 0.323,  "No DBD-specific enrichment",
  53, "GCAGTGAATTATTTACCTCT", 0.098,  "No DBD-specific enrichment",
  54, "CCCAGAAATCATTAACCAGC", 0.107,  "No DBD-specific enrichment",
  55, "GCAGGTTCATTACTAACAGA", 1.940,  "Tier 2 (reproducible DBD-specific)",
  56, "AACATTAAACATTAAACAGT", 0.313,  "No DBD-specific enrichment",
  57, "CCTAGTTAATAATTTGCATC", 1.687,  "Tier 2 (reproducible DBD-specific)",
  58, "TGAGTTAATATTTAGCCCAG", -0.065, "Self-activator (enriched both conditions)",
  59, "TTTTTTAATGGTTAGCCTTT", 0.689,  "No DBD-specific enrichment",
  60, "CATTGTTAATAATTAATACT", 0.985,  "No DBD-specific enrichment",
  61, "CCGATTAACCATTAACCCCC", -0.329, "No DBD-specific enrichment",
  62, "AAAAGCTAATAATTGACAGC", -1.037, "Self-activator (enriched both conditions)",
  63, "GTCTATTAATAATTGACAAA", -1.936, "Self-activator (enriched both conditions)"
)

cat("Reem sequences loaded:", nrow(reem_df), "\n")

# ── HELPER: extract hits from a DNAStringSet ──────────────────────────────────
extract_hits <- function(fa, query, max_mm = 0, strand_label) {
  query_dna <- DNAString(query)
  hits_fwd  <- vmatchPattern(query_dna, fa, max.mismatch = max_mm, fixed = FALSE)
  hits_rev  <- vmatchPattern(reverseComplement(query_dna), fa, max.mismatch = max_mm, fixed = FALSE)
  result_list <- list()
  
  for (strand_hits in list(list(hits_fwd, "+"), list(hits_rev, "-"))) {
    h <- strand_hits[[1]]; s <- strand_hits[[2]]
    for (i in seq_along(h)) {
      m <- h[[i]]
      if (length(m) == 0) next
      result_list[[length(result_list) + 1]] <- data.frame(
        peak_id = names(fa)[i],
        strand  = s,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(result_list) == 0) return(NULL)
  do.call(rbind, result_list)
}

# ── SEARCH FUNCTION ───────────────────────────────────────────────────────────
search_all_reem <- function(fa_path, tissue_label) {
  cat("Loading FASTA:", fa_path, "\n")
  fa <- readDNAStringSet(fa_path)
  cat("  Peaks loaded:", length(fa), "\n")
  
  all_hits <- list()
  for (i in seq_len(nrow(reem_df))) {
    row <- reem_df[i, ]
    h <- extract_hits(fa, row$sequence, max_mm = 0, strand_label = "+")
    if (!is.null(h)) {
      h$reem_id  <- row$reem_id
      h$reem_seq <- row$sequence
      h$delta    <- row$delta
      h$tier     <- row$tier
      h$tissue   <- tissue_label
      all_hits[[length(all_hits) + 1]] <- h
    }
  }
  
  if (length(all_hits) == 0) {
    cat("  No exact matches found for", tissue_label, "\n")
    return(NULL)
  }
  
  result <- do.call(rbind, all_hits)
  cat("  Exact matches:", nrow(result), "| Unique Reem seqs:", length(unique(result$reem_id)), "\n")
  result
}

# ── RUN SEARCHES ──────────────────────────────────────────────────────────────
islet_hits  <- search_all_reem(islet_fa_path,  "Islet")
hepg2_hits  <- search_all_reem(hepg2_fa_path,  "HepG2")

all_hits <- rbind(islet_hits, hepg2_hits)
cat("\nTotal hits (both tissues):", nrow(all_hits), "\n")

# ── PULL PEAK SCORES ──────────────────────────────────────────────────────────
# The CSV uses "chr1:12345-12345" format for peak_id; FASTA uses "chr1_12345-12345"
# We need to harmonise: replace "_" with ":" for the first separator only
harmonise_peak_id <- function(x) {
  # "chr1_12345-12345" → "chr1:12345-12345"
  sub("_", ":", x)
}

islet_csv  <- read_csv(islet_csv_path,  show_col_types = FALSE)
hepg2_csv  <- read_csv(hepg2_csv_path,  show_col_types = FALSE)

peak_scores <- bind_rows(
  islet_csv  %>% select(peak_id, peak_score) %>% mutate(tissue = "Islet"),
  hepg2_csv  %>% select(peak_id, peak_score) %>% mutate(tissue = "HepG2")
)

# Harmonise the hit peak IDs for joining
all_hits$peak_id_csv <- harmonise_peak_id(all_hits$peak_id)

results <- all_hits %>%
  left_join(peak_scores, by = c("peak_id_csv" = "peak_id", "tissue" = "tissue"))

cat("\n=== FINAL RESULTS ===\n")
print(results %>%
        select(tissue, reem_id, reem_seq, peak_id, peak_score, delta, tier) %>%
        arrange(tissue, delta))

# Save
write_csv(results, "/Users/Manu/Downloads/GSE206240_RAW/HNF1A/reem_peak_matches.csv")
cat("\nSaved to: reem_peak_matches.csv\n")

# ── PLOT ──────────────────────────────────────────────────────────────────────
library(ggplot2)

# Tier colour scale (ordered)
tier_colours <- c(
  "Tier 1 (bulletproof DBD-specific)"   = "#d73027",
  "Tier 2 (reproducible DBD-specific)"  = "#fc8d59",
  "No DBD-specific enrichment"           = "#91bfdb",
  "Self-activator (enriched both conditions)" = "#969696"
)

p <- ggplot(results, aes(x = peak_score, y = delta, colour = tier, shape = tissue)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = tier_colours, name = "Reem B1H tier") +
  scale_shape_manual(values = c(Islet = 16, HepG2 = 17), name = "Tissue") +
  labs(
    x     = "ChIP-seq peak score (from BED file)",
    y     = expression(Delta~"DBD – Lib (B1H enrichment score)"),
    title = "HNF1A B1H enrichment vs. ChIP-seq peak score",
    subtitle = "Exact 20bp sequence matches between Reem's B1H sequences and ChIP-seq peaks"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

print(p)

ggsave("/Users/Manu/Downloads/GSE206240_RAW/HNF1A/reem_peak_scatter.pdf",
       plot = p, width = 7, height = 5)
cat("Plot saved to: reem_peak_scatter.pdf\n")