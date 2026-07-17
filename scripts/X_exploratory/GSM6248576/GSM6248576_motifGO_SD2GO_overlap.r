# this script is to compare my motif GO analysis (which uses protein-coding genes) against the paper GO

library(dplyr)
library(ggplot2)
library(VennDiagram)
library(grid)

# convert motif GO result to dataframe
motif_GO <- as.data.frame(GSM6248576_motif)

# paper GO dataframe
paper_GO <- GSM6248576_SD2_GO

# clean paper GO IDs
paper_ids <- unique(na.omit(trimws(as.character(paper_GO$ID))))

# add overlap label + plotting column
motif_GO <- motif_GO %>%
  mutate(
    in_paper = ID %in% paper_ids,
    neglog10_padj = -log10(p.adjust + 1e-300))

# overlapping terms only
motif_GO_overlap <- motif_GO %>%
  filter(in_paper) %>%
  arrange(p.adjust)

# unique terms only
motif_GO_unique <- motif_GO %>%
  filter(!in_paper) %>%
  arrange(p.adjust)

# top 30 terms for plotting
plot_GO_motif <- motif_GO %>%
  arrange(p.adjust) %>%
  slice_head(n = 100) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

# bar plot
ggplot(plot_GO_motif, aes(x = neglog10_padj, y = Description, fill = in_paper)) +
  geom_col() +
  theme_minimal() +
  xlab("-log10 adjusted p-value") +
  ylab("GO Biological Process") +
  labs(fill = "Also In Paper")

# venn diagram
# extract GO IDs
motif_ids <- unique(na.omit(motif_GO$ID))
paper_ids <- unique(na.omit(trimws(as.character(paper_GO$ID))))

# shared / unique
motif_only <- setdiff(motif_ids, paper_ids)
paper_only <- setdiff(paper_ids, motif_ids)
shared_ids <- intersect(motif_ids, paper_ids)

# venn diagram
venn.plot <- venn.diagram(
  x = list(
    Motif_GO = motif_ids,
    Paper_GO = paper_ids),
  category.names = c("Motif + Nearest Protein-Coding GO", "Paper GO"),
  filename = NULL,
  fill = c("skyblue", "pink"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  margin = 0.1)

grid.newpage()
grid.draw(venn.plot)