# performing GO analysis on different gene sets derived from peaks
# what biological processes are enriched among these peak-associated genes?

# load libraries
library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# promoter-proximal genes (TSS' found within the 10kb window of the peak)
GSM6248576_promoterproximal_genes <- unique(unlist(strsplit(GSM6248576$genes_TSS_10kb[GSM6248576$genes_TSS_10kb != ""],",")))

# gene-body overlap genes
GSM6248576_body_genes <- unique(unlist(strsplit(GSM6248576$genes_geneBodyOverlap_10kb[GSM6248576$genes_geneBodyOverlap_10kb != ""],",")))

# nearest genes
GSM6248576_nearest_genes <- unique(GSM6248576$nearest_gene)

# nearest protein coding genes
GSM6248576_nearest_proteincoding_genes <- unique(trimws(GSM6248576$nearest_proteincoding_gene))

# promoter genes (genes with promoter annotation)
GSM6248576_promoter_genes <- GSM6248576 %>%
  filter(!is.na(peak_overlap_type),
         grepl("promoter", peak_overlap_type, ignore.case = TRUE)) %>%
  pull(nearest_gene) %>%
  trimws() %>%
  .[!is.na(.) & . != ""] %>%
  unique()

# nearest protein coding + promoter genes
GSM6248576_proteincoding_promoter_genes <- GSM6248576 %>%
  filter(str_detect(peak_overlap_type, "Promoter"),
         !is.na(nearest_proteincoding_gene)) %>%
  pull(nearest_proteincoding_gene) %>%
  unique()


# FIMO genes (genes that get a motif hit on FIMO)

GSM6248576_FIMO_peaks_fixed <- GSM6248576_FIMO$sequence_name %>%
  unique() %>%
  gsub("_", ":", .)

# subset main peak table to only motif-containing peaks
GSM6248576_motif_peaks <- GSM6248576 %>%
  filter(peak_id %in% GSM6248576_FIMO_peaks_fixed)

# extract nearest protein-coding genes from those motif-containing peaks
GSM6248576_motif_genes <- GSM6248576_motif_peaks %>%
  pull(nearest_proteincoding_gene) %>%
  trimws() %>%
  .[!is.na(.) & . != ""] %>%
  unique()

# helper function to convert gene SYMBOLs to ENTREZ IDs
symbol_to_entrez <- function(symbols) {
  symbols <- unique(trimws(symbols))
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  bitr(symbols,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db) |>
    dplyr::pull(ENTREZID) |>
    unique()}

# apply to all gene sets
GSM6248576_promoterproximal_entrez      <- symbol_to_entrez(GSM6248576_promoterproximal_genes)
GSM6248576_body_entrez                  <- symbol_to_entrez(GSM6248576_body_genes)
GSM6248576_nearest_entrez               <- symbol_to_entrez(GSM6248576_nearest_genes)
GSM6248576_nearest_proteincoding_entrez <- symbol_to_entrez(GSM6248576_nearest_proteincoding_genes)
GSM6248576_promoter_entrez              <- symbol_to_entrez(GSM6248576_promoter_genes)
GSM6248576_motif_entrez                 <- symbol_to_entrez(GSM6248576_motif_genes)
GSM6248576_proteincoding_promoter_entrez <- symbol_to_entrez(GSM6248576_proteincoding_promoter_genes)

# run GO enrichment on proximal-promoter genes
GSM6248576_promoterproximal <- enrichGO(
  gene          = GSM6248576_promoterproximal_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

View(as.data.frame(GSM6248576_promoterproximal))

# run GO enrichment on genebody genes
GSM6248576_body_pvalue <- enrichGO(
  gene          = GSM6248576_body_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

View(as.data.frame(GSM6248576_body_pvalue))

# run GO enrichment on nearest genes
GSM6248576_nearest <- enrichGO(
  gene          = GSM6248576_nearest_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

View(as.data.frame(GSM6248576_nearest))

# run GO enrichment on nearest protein-coding genes
GSM6248576_nearest_proteincoding <- enrichGO(
  gene          = GSM6248576_nearest_proteincoding_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

View(as.data.frame(GSM6248576_nearest_proteincoding))

# run GO enrichment on promoter genes
GSM6248576_promoter <- enrichGO(
  gene          = GSM6248576_promoter_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

View(as.data.frame(GSM6248576_promoter))

GSM6248576_promoter_all <- enrichGO(
  gene          = GSM6248576_promoter_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  readable      = TRUE)

View(as.data.frame(GSM6248576_promoter_all))

# run GO enrichment on protein-coding + promoter genes
GSM6248576_GO_proteincoding_promoter <- enrichGO(
  gene          = GSM6248576_proteincoding_promoter_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

View(as.data.frame(GSM6248576_GO_proteincoding_promoter))

# run GO enrichment on motif genes
GSM6248576_motif <- enrichGO(
  gene          = GSM6248576_motif_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

View(as.data.frame(GSM6248576_motif))
