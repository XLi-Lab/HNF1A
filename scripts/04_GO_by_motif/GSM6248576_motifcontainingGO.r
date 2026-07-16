# this r script computes GO BP enrichment for the nearest protein-coding genes of
# motif-containing islet peaks (the motif-restricted gene set for the Figure 7 comparison)
# output saved as GSM6248576_ProteinCodingGO_Motif.csv

# load libraries
library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# read the per-peak table (from motif scanning)
all_peaks <- read_csv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_FIMO_AllPeaks.csv", show_col_types = FALSE)

# helper: gene SYMBOLs -> ENTREZ IDs
symbol_to_entrez <- function(symbols) {
  symbols <- unique(trimws(symbols))
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) |>
    dplyr::pull(ENTREZID) |>
    unique()}

# motif-containing peaks -> nearest protein-coding genes -> ENTREZ
motif_genes  <- all_peaks %>%
  filter(motif_group == "Motif") %>%
  pull(nearest_proteincoding_gene) %>%
  trimws() %>%
  .[!is.na(.) & . != ""] %>%
  unique()

motif_entrez <- symbol_to_entrez(motif_genes)

# GO BP enrichment on motif-containing gene set
GSM6248576_GO_motif <- enrichGO(
  gene          = motif_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

# write the motif-containing GO for the Figure 7 comparison
write.csv(as.data.frame(GSM6248576_GO_motif), "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_ProteinCodingGO_Motif.csv", row.names = FALSE)