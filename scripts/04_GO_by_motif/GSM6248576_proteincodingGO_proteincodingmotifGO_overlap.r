# this r script is to compare the protein coding gene list GO analysis compared to the protein coding + motif GO analysis

# load libraries
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(grid)

# convert enrichGO results to dataframes
go_pc       <- read.csv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_ProteinCodingGO_AllPeaks.csv")
go_pc_motif <- read.csv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_ProteinCodingGO_Motif.csv")

# get GO IDs
pc_ids <- unique(na.omit(trimws(go_pc$ID)))
motif_ids <- unique(na.omit(trimws(go_pc_motif$ID)))

# label overlaps between the 2 analyses
go_pc <- go_pc %>%
  mutate(
    in_motif = ID %in% motif_ids,
    neglog10_padj = -log10(p.adjust))

go_pc_motif <- go_pc_motif %>%
  mutate(
    in_pc = ID %in% pc_ids,
    neglog10_padj = -log10(p.adjust))

# counts
sum(go_pc$in_motif)         # overlap
sum(!go_pc$in_motif)        # only in protein-coding

sum(go_pc_motif$in_pc)      # overlap
sum(!go_pc_motif$in_pc)     # only in protein-coding + motif

# terms lost after motif filtering
pc_only_ids <- setdiff(pc_ids, motif_ids)
pc_only_terms <- go_pc %>%
  filter(ID %in% pc_only_ids) %>%
  arrange(p.adjust)

# terms unique to motif subset
motif_only_ids <- setdiff(motif_ids, pc_ids)
motif_only_terms <- go_pc_motif %>%
  filter(ID %in% motif_only_ids) %>%
  arrange(p.adjust)

# overlapping terms
overlap_ids <- intersect(pc_ids, motif_ids)

pc_overlap_terms <- go_pc %>%
  filter(ID %in% overlap_ids) %>%
  arrange(p.adjust)

motif_overlap_terms <- go_pc_motif %>%
  filter(ID %in% overlap_ids) %>%
  arrange(p.adjust)

# top 30 protein-coding terms for plotting
plot_pc <- go_pc %>%
  arrange(p.adjust) %>%
  slice_head(n = 30) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

ggplot(plot_pc, aes(x = neglog10_padj, y = Description, fill = in_motif)) +
  geom_col() +
  theme_minimal() +
  xlab("-log10 adjusted p-value") +
  ylab("GO Biological Process") +
  labs(fill = "Also In Motif Subset")

# top 30 motif terms for plotting
plot_motif <- go_pc_motif %>%
  arrange(p.adjust) %>%
  slice_head(n = 30) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

ggplot(plot_motif, aes(x = neglog10_padj, y = Description, fill = in_pc)) +
  geom_col() +
  theme_minimal() +
  xlab("-log10 adjusted p-value") +
  ylab("GO Biological Process") +
  labs(fill = "Also In Protein-Coding Set")

# venn diagram
venn.plot <- venn.diagram(
  x = list(
    protein_coding = pc_ids,
    protein_coding_motif = motif_ids),
  filename = NULL,
  disable.logging = TRUE,
  fill = c("skyblue", "pink"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  main = "GO ID overlap")

grid.newpage()
grid.draw(venn.plot)

