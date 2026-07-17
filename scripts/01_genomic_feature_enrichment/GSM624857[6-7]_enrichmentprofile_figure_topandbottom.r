# this r script builds a two-panel figure (genomic-feature fold enrichment) for islet and HepG2 HNF1A peaks

# load libraries
library(ggplot2)
library(dplyr)
library(patchwork)

# read the per-tissue enrichment results written by the profile scripts
islet <- read.csv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_Enrichment_Profile.csv")
hep   <- read.csv("HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_Enrichment_Profile.csv")

# relabel UTR categories with prime notation for display
relabel_utr <- function(x) {
  x[x == "5UTR"] <- "5' UTR"
  x[x == "3UTR"] <- "3' UTR" 
  x}
islet$category <- relabel_utr(islet$category)
hep$category   <- relabel_utr(hep$category)

obs_islets      <- setNames(islet$obs,      islet$category)
exp_prob_islets <- setNames(islet$exp_prob, islet$category)
obs_hepg2       <- setNames(hep$obs,        hep$category)
exp_prob_hepg2  <- setNames(hep$exp_prob,   hep$category)

# build fe_df for islets
fe_islets <- data.frame(
  category       = names(obs_islets),
  obs_n          = as.integer(obs_islets),
  total_n        = sum(obs_islets),
  FoldEnrichment = (as.numeric(obs_islets) / sum(obs_islets)) / exp_prob_islets,
  cell_type      = "Pancreatic Islets  (n = 426)")

# build fe_df for HepG2
fe_hepg2 <- data.frame(
  category       = names(obs_hepg2),
  obs_n          = as.integer(obs_hepg2),
  total_n        = sum(obs_hepg2),
  FoldEnrichment = (as.numeric(obs_hepg2) / sum(obs_hepg2)) / exp_prob_hepg2,
  cell_type      = "HepG2  (n = 4,784)")

# combine and add percentage label
fe_combined <- rbind(fe_islets, fe_hepg2) |>
  mutate(
    pct       = round(100 * obs_n / total_n, 1),
    bar_label = paste0("n=", obs_n, " (", pct, "%)"),
    cell_type = factor(cell_type,
                       levels = c("Pancreatic Islets  (n = 426)",
                                  "HepG2  (n = 4,784)")))

# enforce consistent category ordering (by islet fold enrichment)
cat_order <- fe_islets |>
  arrange(FoldEnrichment) |>
  pull(category)

fe_combined$category <- factor(fe_combined$category, levels = cat_order)

# panel colours
cell_colours <- c(
  "Pancreatic Islets  (n = 426)" = "#1a80bb",
  "HepG2  (n = 4,784)"           = "#a00000")

# shared theme
enrichment_theme <- theme_grey(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.y        = element_text(size = 11, face = "bold"),
    axis.title.x       = element_text(face = "bold", size = 12),
    axis.text.x        = element_text(face = "bold", size = 11))

# islets plot
p_islets <- ggplot(fe_combined |> filter(cell_type == "Pancreatic Islets  (n = 426)"),
                   aes(x = category, y = FoldEnrichment, fill = cell_type)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#B8B8B8") +
  geom_text(aes(label = paste0("n=", obs_n, "\n(", pct, "%)")),
            hjust = -0.05, size = 3.6, lineheight = 0.85, colour = "grey25", fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = cell_colours, guide = "none") +
  scale_y_continuous(
    limits = c(0, 8.25),
    breaks = seq(0, 8, by = 1),
    expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Pancreatic Islet Peaks",
    x     = NULL,
    y     = "Fold Enrichment (Observed/Expected)") +
  enrichment_theme

# HepG2 plot
p_hepg2 <- ggplot(fe_combined |> filter(cell_type == "HepG2  (n = 4,784)"),
                  aes(x = category, y = FoldEnrichment, fill = cell_type)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#B8B8B8") +
  geom_text(aes(label = paste0("n=", obs_n, "\n(", pct, "%)")),
            hjust = -0.08, size = 3.6, lineheight = 0.85, colour = "grey25", fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = cell_colours, guide = "none") +
  scale_y_continuous(
    limits = c(0, 8.25),
    breaks = seq(0, 8, by = 1),
    expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "HepG2 Peaks",
    x     = NULL,
    y     = "Fold Enrichment (Observed/Expected)") +
  enrichment_theme

(p_islets / p_hepg2) +
  plot_annotation(
    title = "Genomic Feature Enrichment of HNF1A ChIP-Seq Peaks",
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5)))

ggsave("enrichment_stacked.pdf", width = 8, height = 12, dpi = 300)
