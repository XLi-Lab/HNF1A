# Extract the 14-bp core from Reem's sequence (positions 2-15)
reem_core <- "GGTTAATAATTAAC"

# Search in the combined FIMO object
matches <- GSM6248576.7_Combined_FIMO %>%
  filter(matched_sequence == reem_core)

table(matches$tissue)



# search the full 17bp

library(Biostrings)

# Load the full peak sequences
islet_fa  <- readDNAStringSet("/Users/Manu/Downloads/GSE206240_RAW/HNF1A/HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_peaks-FASTAs.fa")
hepg2_fa  <- readDNAStringSet("/Users/Manu/Downloads/GSE206240_RAW/HNF1A/HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_peaks-FASTAs.fa")

query <- "TTGGTTAATAATTAACT"

islet_hits <- vmatchPattern(query, islet_fa, fixed = TRUE)
hep_hits   <- vmatchPattern(query, hepg2_fa, fixed = TRUE)

# Also check reverse complement
query_rc <- as.character(reverseComplement(DNAString(query)))
islet_hits_rc <- vmatchPattern(query_rc, islet_fa, fixed = TRUE)
hep_hits_rc   <- vmatchPattern(query_rc, hepg2_fa, fixed = TRUE)

cat("Islet fwd:", sum(lengths(islet_hits)), "RC:", sum(lengths(islet_hits_rc)), "\n")
cat("HepG2 fwd:", sum(lengths(hep_hits)),   "RC:", sum(lengths(hep_hits_rc)),   "\n")


# search the full 17bp with 1 variaton

library(Biostrings)

# ── paths ──────────────────────────────────────────────────────────────────────
islet_fa_path <- "/Users/Manu/Downloads/GSE206240_RAW/HNF1A/HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_peaks-FASTAs.fa"
hepg2_fa_path <- "/Users/Manu/Downloads/GSE206240_RAW/HNF1A/HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_peaks-FASTAs.fa"

islet_fa <- readDNAStringSet(islet_fa_path)
hepg2_fa <- readDNAStringSet(hepg2_fa_path)

query    <- DNAString("TTGGTTAATAATTAACT")
query_rc <- reverseComplement(query)

# ── robust hit extractor ───────────────────────────────────────────────────────
extract_hits <- function(fa, q, max_mm = 1, strand_label) {
  hits <- vmatchPattern(q, fa, max.mismatch = max_mm, fixed = FALSE)
  result_list <- list()
  for (i in seq_along(hits)) {
    m <- hits[[i]]
    if (length(m) == 0) next
    for (j in seq_along(m)) {
      s <- start(m)[j]
      e <- end(m)[j]
      matched_seq <- as.character(subseq(fa[[i]], start = s, end = e))
      result_list[[length(result_list) + 1]] <- data.frame(
        peak    = names(fa)[i],
        start   = s,
        end     = e,
        strand  = strand_label,
        matched = matched_seq,
        stringsAsFactors = FALSE)
    }
  }
  if (length(result_list) == 0) {
    return(data.frame(peak = character(), start = integer(),
                      end = integer(), strand = character(),
                      matched = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, result_list)
}

# ── run ────────────────────────────────────────────────────────────────────────
islet_fwd <- extract_hits(islet_fa, query,    max_mm = 1, strand_label = "+")
islet_rev <- extract_hits(islet_fa, query_rc, max_mm = 1, strand_label = "-")
hepg2_fwd <- extract_hits(hepg2_fa, query,    max_mm = 1, strand_label = "+")
hepg2_rev <- extract_hits(hepg2_fa, query_rc, max_mm = 1, strand_label = "-")

islet_hits <- rbind(islet_fwd, islet_rev)
hepg2_hits <- rbind(hepg2_fwd, hepg2_rev)

# ── summarise ─────────────────────────────────────────────────────────────────
cat("=== ISLETS (GSM6248576) ===\n")
cat("Total hits (<=1 mismatch):", nrow(islet_hits), "\n")
if (nrow(islet_hits) > 0) print(islet_hits)

cat("\n=== HepG2 (GSM6248577) ===\n")
cat("Total hits (<=1 mismatch):", nrow(hepg2_hits), "\n")
if (nrow(hepg2_hits) > 0) print(hepg2_hits)