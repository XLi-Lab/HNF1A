# this r script is to make the HNF1A motif-containg vs motif-lacking proportion plots,
# for both cell types: primary islets (GSM6248576) and hepG2 (GSM6248577)

# load libraries
library(dplyr)
library(ggplot2)

combined <- bind_rows(
  GSM6248576_FIMO_AllPeaks %>% mutate(cell_type = "Pancreatic Islets"),
  GSM6248577_FIMO_AllPeaks %>% mutate(cell_type = "HepG2")) %>%
  mutate(cell_type = factor(cell_type,
                            levels = c("Pancreatic Islets", "HepG2")))

summary_tbl <- combined %>%
  count(cell_type, motif_group) %>%
  group_by(cell_type) %>%
  mutate(
    total      = sum(n),
    proportion = n / total,
    pct_label  = paste0(round(proportion * 100, 1), "%\n(n = ", n, ")")) %>%
  ungroup() %>%
  mutate(
    fill_colour = case_when(
      cell_type == "Pancreatic Islets" & motif_group == "Motif"    ~ "#1a80bb",
      cell_type == "Pancreatic Islets" & motif_group == "No Motif" ~ "#c9cdd6",
      cell_type == "HepG2"             & motif_group == "Motif"    ~ "#a00000",
      cell_type == "HepG2"             & motif_group == "No Motif" ~ "#c9cdd6"))

p <- ggplot(summary_tbl,
            aes(x = "", y = proportion, fill = fill_colour)) +
  geom_col(width = 1, colour = "white", linewidth = 0.5) +
  geom_text(aes(label = pct_label),
          position = position_stack(vjust = 0.5),
          size = 6, fontface = "bold", lineheight = 1.1,
          colour = ifelse(summary_tbl$motif_group == "Motif", "white", "#444444")) +
  coord_polar(theta = "y") +
  scale_fill_identity(
    guide  = "legend",
    labels = c("#1a80bb" = "With HNF1A Motif (Islets)",
               "#a00000" = "With HNF1A Motif (HepG2)",
               "#c9cdd6" = "No HNF1A Motif"),
    breaks = c("#1a80bb", "#a00000", "#c9cdd6")) +
  facet_wrap(~ cell_type, labeller = as_labeller(c(
    "Pancreatic Islets" = "Islets\n(n = 426)",
    "HepG2"             = "HepG2\n(n = 4,784)"))) +
  labs(
    title = "Proportion of HNF1A ChIP-Seq Peaks Containing A Canonical HNF1A Motif",
    fill  = NULL) +
  theme_void(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 18,
                                   margin = margin(b = 10, t = 10)),
    strip.text = element_text(face = "bold", size = 18, margin = margin(t = 20)),
    legend.position = "bottom",
    legend.text     = element_text(size = 16, face = "bold"),
    legend.key.spacing = unit(0.5, "cm"))

print(p)

ggsave("HNF1A_motif_proportion_pie.pdf", plot = p, width = 7, height = 5, device = "pdf")
ggsave("HNF1A_motif_proportion_pie.png", plot = p, width = 7, height = 5, dpi = 300)