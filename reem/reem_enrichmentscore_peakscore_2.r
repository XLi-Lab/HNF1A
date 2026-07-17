# ============================================================
# Reem × ChIP-seq sequence matching  (v2)
# 
# BEFORE RUNNING:
#   1. Open HNF1A_B1H_all_motifs_corrected_2.numbers in Numbers
#   2. File → Export To → CSV → save as "reem_b1h.csv"
#      in the same folder as this script, or update reem_csv_path below
#
# The script does:
#   (a) Exact 20bp match
#   (b) 14bp core match (shave 3bp flanking from each end)
#   Both against islet + HepG2 peak FASTAs.
#   Peak scores are pulled from the FIMO_AllPeaks CSVs.
#   Two scatter plots are produced (one per match strategy).
# ============================================================

library(Biostrings)
library(readr)
library(dplyr)
library(ggplot2)

# ── PATHS — edit if needed ────────────────────────────────────────────────────
base_dir <- "/Users/Manu/Downloads/GSE206240_RAW/HNF1A"

reem_csv_path <- "/Users/Manu/Desktop/HNF1A_B1H_all_motifs_corrected_2.csv"

islet_fa_path  <- file.path(base_dir, "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_peaks-FASTAs.fa")
islet_csv_path <- file.path(base_dir, "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_FIMO_AllPeaks.csv")

hepg2_fa_path  <- file.path(base_dir, "HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_peaks-FASTAs.fa")
hepg2_csv_path <- file.path(base_dir, "HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_FIMO_AllPeaks.csv")

# ── LOAD REEM'S DATA ──────────────────────────────────────────────────────────
reem_raw <- read_csv(reem_csv_path, show_col_types = FALSE)

# Expect columns: motif_id, sequence, delta_DBD_minus_Lib, tier  (among others)
# Rename for convenience
reem_df <- reem_raw %>%
  select(
    reem_id  = motif_id,
    sequence = sequence,
    delta    = delta_DBD_minus_Lib,
    tier     = tier
  ) %>%
  mutate(
    sequence = toupper(trimws(sequence)),
    reem_id  = as.integer(reem_id),
    # Derive 14bp core: remove 3bp from each end
    # Only valid for 20bp sequences; the one 17bp seq gets 3bp trimmed each side → 11bp
    core = substr(sequence, 4, nchar(sequence) - 3)
  )

cat("Reem sequences loaded:", nrow(reem_df), "\n")
cat("Sequence lengths present:", paste(unique(nchar(reem_df$sequence)), collapse = ", "), "\n")
cat("Core lengths present:",     paste(unique(nchar(reem_df$core)),     collapse = ", "), "\n\n")

# ── LOAD PEAK SCORES ──────────────────────────────────────────────────────────
# The FIMO CSVs use "chr1:start-end" format for peak_id;
# the FASTA headers use "chr1_start-end". We harmonise by replacing the first "_"
harmonise_id <- function(x) sub("_", ":", x)

islet_scores <- read_csv(islet_csv_path, show_col_types = FALSE) %>%
  select(peak_id, peak_score) %>%
  mutate(tissue = "Islet")

hepg2_scores <- read_csv(hepg2_csv_path, show_col_types = FALSE) %>%
  select(peak_id, peak_score) %>%
  mutate(tissue = "HepG2")

peak_scores <- bind_rows(islet_scores, hepg2_scores)

# ── SEARCH HELPER ─────────────────────────────────────────────────────────────
# Returns a data frame of (peak_id, strand) for all peaks containing `query`
# (or its reverse complement) with at most `max_mm` mismatches (0 = exact).
search_fa <- function(fa, query, max_mm = 0) {
  q     <- DNAString(query)
  q_rc  <- reverseComplement(q)
  
  hits_fwd <- vmatchPattern(q,    fa, max.mismatch = max_mm, fixed = FALSE)
  hits_rev <- vmatchPattern(q_rc, fa, max.mismatch = max_mm, fixed = FALSE)
  
  collect <- function(hits, strand_label) {
    idx <- which(lengths(hits) > 0)
    if (length(idx) == 0) return(NULL)
    data.frame(peak_id = names(fa)[idx], strand = strand_label,
               stringsAsFactors = FALSE)
  }
  
  rbind(collect(hits_fwd, "+"), collect(hits_rev, "-"))
}

# Search all Reem sequences against one FASTA; `query_col` is "sequence" or "core"
search_all <- function(fa_path, tissue_label, reem_df, query_col) {
  cat("Loading", tissue_label, "FASTA ...\n")
  fa <- readDNAStringSet(fa_path)
  cat("  Peaks:", length(fa), "\n")
  
  results <- lapply(seq_len(nrow(reem_df)), function(i) {
    row   <- reem_df[i, ]
    query <- row[[query_col]]
    h     <- search_fa(fa, query, max_mm = 0)
    if (is.null(h)) return(NULL)
    h$reem_id  <- row$reem_id
    h$reem_seq <- row$sequence
    h$query    <- query
    h$delta    <- row$delta
    h$tier     <- row$tier
    h$tissue   <- tissue_label
    h
  })
  
  out <- do.call(rbind, Filter(Negate(is.null), results))
  n_hits <- if (is.null(out)) 0 else nrow(out)
  n_reem <- if (is.null(out)) 0 else length(unique(out$reem_id))
  cat("  Hits:", n_hits, "| Unique Reem seqs matched:", n_reem, "\n")
  out
}

# ── JOIN PEAK SCORES ──────────────────────────────────────────────────────────
add_peak_scores <- function(hits_df) {
  if (is.null(hits_df) || nrow(hits_df) == 0) return(hits_df)
  hits_df %>%
    mutate(peak_id_csv = harmonise_id(peak_id)) %>%
    left_join(peak_scores, by = c("peak_id_csv" = "peak_id", "tissue" = "tissue")) %>%
    select(tissue, reem_id, reem_seq, query, peak_id, peak_score, delta, tier, strand)
}

# ── RUN: EXACT 20bp ───────────────────────────────────────────────────────────
cat("\n=== EXACT 20bp SEARCH ===\n")
exact_islet <- search_all(islet_fa_path, "Islet", reem_df, "sequence")
exact_hepg2 <- search_all(hepg2_fa_path, "HepG2", reem_df, "sequence")
exact_all   <- add_peak_scores(rbind(exact_islet, exact_hepg2))

cat("\nExact 20bp results:\n")
print(exact_all %>% arrange(tissue, delta) %>%
        select(tissue, reem_id, delta, tier, peak_id, peak_score))

# ── RUN: 14bp CORE ────────────────────────────────────────────────────────────
cat("\n=== 14bp CORE SEARCH ===\n")
core_islet <- search_all(islet_fa_path, "Islet", reem_df, "core")
core_hepg2 <- search_all(hepg2_fa_path, "HepG2", reem_df, "core")
core_all   <- add_peak_scores(rbind(core_islet, core_hepg2))

cat("\n14bp core results:\n")
print(core_all %>% arrange(tissue, delta) %>%
        select(tissue, reem_id, delta, tier, peak_id, peak_score, query))

# ── SAVE CSVs ─────────────────────────────────────────────────────────────────
write_csv(exact_all, file.path(out_dir, "reem_exact20bp_matches.csv"))
write_csv(core_all,  file.path(out_dir, "reem_core14bp_matches.csv"))
cat("\nCSVs saved.\n")

# ── PLOT HELPER ───────────────────────────────────────────────────────────────
tier_colours <- c(
  "Tier 1 (bulletproof DBD-specific)"        = "#d73027",
  "Tier 2 (reproducible DBD-specific)"       = "#fc8d59",
  "No DBD-specific enrichment"               = "#91bfdb",
  "Self-activator (enriched both conditions)"= "#969696"
)

make_plot <- function(df, title_suffix) {
  if (is.null(df) || nrow(df) == 0) {
    message("No data for plot: ", title_suffix); return(NULL)
  }
  ggplot(df, aes(x = peak_score, y = delta, colour = tier, shape = tissue)) +
    geom_point(size = 3, alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = tier_colours, name = "B1H tier") +
    scale_shape_manual(values = c(Islet = 16, HepG2 = 17), name = "Tissue") +
    labs(
      x        = "ChIP-seq peak score",
      y        = expression(Delta ~ "DBD \u2013 Lib (B1H enrichment score)"),
      title    = paste("HNF1A B1H enrichment vs. ChIP-seq peak score"),
      subtitle = title_suffix
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
}

p_exact <- make_plot(exact_all, "Exact 20bp match")
p_core  <- make_plot(core_all,  "14bp core match (±3bp flanking removed)")

if (!is.null(p_exact)) {
  print(p_exact)
  ggsave(file.path(out_dir, "reem_scatter_exact20bp.pdf"), p_exact, width = 7, height = 5)
}
if (!is.null(p_core)) {
  print(p_core)
  ggsave(file.path(out_dir, "reem_scatter_core14bp.pdf"),  p_core,  width = 7, height = 5)
}

cat("\nDone. Output files saved to:", out_dir, "\n")