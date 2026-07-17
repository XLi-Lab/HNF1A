library(dplyr)
library(GenomicRanges)
library(tidyr)

# ---- 1. build GRanges of all peaks per tissue from peak_id ----
# peak_id format is "chr1:3219340-3219738"
peakid_to_granges <- function(peak_ids) {
  parts <- do.call(rbind, strsplit(peak_ids, "[:-]"))
  GRanges(seqnames = parts[, 1],
          ranges   = IRanges(start = as.integer(parts[, 2]),
                             end   = as.integer(parts[, 3])),
          peak_id  = peak_ids)
}

islet_gr <- peakid_to_granges(GSM6248576$peak_id)
hep_gr   <- peakid_to_granges(GSM6248577$peak_id)

# ---- 2. coordinate-based tissue-specificity ----
overlap_classify <- function(gr_self, gr_other, minoverlap = 1L) {
  hits <- findOverlaps(gr_self, gr_other, minoverlap = minoverlap)
  shared_idx <- unique(queryHits(hits))
  out <- rep("Specific", length(gr_self))
  out[shared_idx] <- "Shared"
  out
}

islet_coord_class <- overlap_classify(islet_gr, hep_gr)
hep_coord_class   <- overlap_classify(hep_gr,   islet_gr)

cat("---- COORDINATE-BASED (any overlap) ----\n")
cat("Islet shared fraction:", round(mean(islet_coord_class == "Shared"), 3), "\n")
cat("HepG2 shared fraction:", round(mean(hep_coord_class   == "Shared"), 3), "\n\n")

# ---- 3. gene-based composition ----
cat("---- GENE-BASED (nearest protein-coding gene) ----\n")
cat("Shared genes:", length(shared_genes),
    "| Islet-only:", length(islet_only),
    "| HepG2-only:", length(hep_only), "\n\n")

# ---- 4. per-peak agreement (islet) ----
islet_compare <- GSM6248576 %>%
  transmute(peak_id,
            gene = trimws(nearest_proteincoding_gene),
            coord_label = islet_coord_class) %>%
  mutate(gene_label = case_when(
    gene %in% shared_genes ~ "Shared",
    gene %in% islet_only   ~ "Specific",
    TRUE                   ~ NA_character_)) %>%
  filter(!is.na(gene_label))

cat("---- AGREEMENT (islet): gene rows vs coord cols ----\n")
print(table(gene = islet_compare$gene_label, coord = islet_compare$coord_label))
cat("\nDisagree:", round(100 * mean(islet_compare$gene_label != islet_compare$coord_label), 1), "%\n")

# ---- 5. overlap-threshold sensitivity ----
cat("\n---- shared fraction vs minimum overlap (islet) ----\n")
for (mo in c(1L, 50L, 100L, 200L)) {
  cl <- overlap_classify(islet_gr, hep_gr, minoverlap = mo)
  cat(sprintf("minoverlap = %4d bp : shared = %.3f\n", mo, mean(cl == "Shared")))
}