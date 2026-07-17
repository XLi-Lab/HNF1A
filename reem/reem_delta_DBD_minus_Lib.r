library(ggplot2)
library(dplyr)
library(tidyr)

# convert both results to dataframes
go_pc    <- as.data.frame(GSM6248576_nearest_proteincoding) %>%
  mutate(neglog10_padj = -log10(p.adjust), subset = "All peaks\n(protein-coding genes)")

go_motif <- as.data.frame(GSM6248576_motif) %>%
  mutate(neglog10_padj = -log10(p.adjust), subset = "Motif-containing peaks\n(direct binding)")

# pick the top 15 terms from the protein-coding GO as your reference set
top_terms <- go_pc %>%
  arrange(p.adjust) %>%
  slice_head(n = 30) %>%
  pull(Description)

# filter both dataframes to only those terms
plot_data <- bind_rows(go_pc, go_motif) %>%
  filter(Description %in% top_terms) %>%
  mutate(
    Description = factor(Description, levels = rev(top_terms)),
    subset = factor(subset, levels = c(
      "All peaks\n(protein-coding genes)",
      "Motif-containing peaks\n(direct binding)"))
  )

# plot
ggplot(plot_data, aes(x = neglog10_padj, y = Description, size = Count, colour = neglog10_padj)) +
  geom_point() +
  facet_wrap(~ subset, ncol = 2) +
  scale_colour_gradient(low = "#56B4E9", high = "#D55E00", name = "-log10\nadj. p-value") +
  scale_size_continuous(name = "Gene count", range = c(2, 8)) +
  theme_bw(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.y      = element_text(size = 9),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "-log10 adjusted p-value",
    y = NULL,
    title = "GO Biological Process enrichment",
    subtitle = "Islet HNF1A peaks — all peaks vs. direct-binding subset"
  )