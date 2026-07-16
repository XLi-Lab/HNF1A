# this r script computes the all-peaks protein-coding GO enrichment and benchmarks it
# against the paper's published GO (arm 2, project-derived gene list)

# load libraries
library(readr)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(VennDiagram)
library(grid)

# read the per-peak table and the published paper's GO
all_peaks <- read_csv("HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_peakcentric.csv", show_col_types = FALSE)
GSM6248577_SD2_GO <- readxl::read_xlsx("HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_SD2_GO.xlsx")

# helper function which converts gene SYMBOLs -> ENTREZ IDs
symbol_to_entrez <- function(symbols) {
  symbols <- unique(trimws(symbols))
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) |>
    dplyr::pull(ENTREZID) |>
    unique()}

# all-peaks nearest protein-coding genes -> ENTREZ -> GO
all_genes  <- all_peaks$nearest_proteincoding_gene
all_entrez <- symbol_to_entrez(all_genes)

GSM6248577_GO_nearest_proteincoding <- enrichGO(
  gene          = all_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

# compare against the paper's published GO
my_GO <- as.data.frame(GSM6248577_GO_nearest_proteincoding)
paper_ids <- unique(na.omit(trimws(GSM6248577_SD2_GO$ID)))

my_GO <- my_GO %>%
  mutate(
    in_paper = ID %in% paper_ids,
    neglog10_padj = -log10(p.adjust))

# counts
sum(my_GO$in_paper)     # overlap with paper
sum(!my_GO$in_paper)    # unique to this analysis

# overlapping terms only
my_GO_overlap <- my_GO %>%
  filter(in_paper) %>%
  arrange(p.adjust)

# unique to this analysis
my_GO_unique <- my_GO %>%
  filter(!in_paper) %>%
  arrange(p.adjust)

# bar plot: top 30 terms, coloured by whether also in paper
plot_GO <- my_GO %>%
  arrange(p.adjust) %>%
  slice_head(n = 30) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

ggplot(plot_GO, aes(x = neglog10_padj, y = Description, fill = in_paper)) +
  geom_col() +
  theme_minimal() +
  xlab("-log10 adjusted p-value") +
  ylab("GO Biological Process") +
  labs(fill = "Also in paper")

# venn diagram to plot concordance between my GO IDs vs paper GO IDs
my_ids <- unique(na.omit(trimws(my_GO$ID)))

grid.newpage()
venn.plot <- venn.diagram(
  x = list(
    `Present Analysis` = my_ids,
    Published = paper_ids),
  filename = NULL,
  fill = c("skyblue", "pink"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "GO ID Overlap",
  disable.logging = TRUE)
grid::grid.draw(venn.plot)

# GO IDs present in paper but absent from this analysis
paper_only_ids <- setdiff(paper_ids, my_ids)
paper_only_terms <- GSM6248577_SD2_GO %>%
  filter(ID %in% paper_only_ids) %>%
  arrange(p.adjust)

# recovery of published terms
recovery_hepatocytes <- length(intersect(my_ids, paper_ids)) / length(paper_ids) * 100
recovery_hepatocytes

# write the all-peaks protein-coding GO
write.csv(my_GO, "HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_ProteinCodingGO_AllPeaks.csv", row.names = FALSE)

