# this r script builds figure 7 (side-by-side GO BP dot plot comparing all islet peaks
# vs motif-containing peaks), highlighting immune/T-cell terms lost on motif restriction

# load libraries
library(ggplot2)
library(dplyr)

go_pc <- read.csv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_ProteinCodingGO_AllPeaks.csv") %>%
  mutate(neglog10_padj = -log10(p.adjust),
         subset = "All Peaks")

go_motif <- read.csv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_ProteinCodingGO_Motif.csv") %>%
  mutate(neglog10_padj = -log10(p.adjust),
         subset = "Motif-Containing\nPeaks")

top_terms <- go_pc %>% filter(p.adjust < 0.02) %>% arrange(p.adjust) %>% pull(Description)

plot_data <- bind_rows(go_pc, go_motif) %>%
  filter(Description %in% top_terms) %>%
  mutate(
    Description = factor(Description, levels = rev(top_terms)),
    subset = factor(subset, levels = c(
      "All Peaks",
      "Motif-Containing\nPeaks")), x = 1)   

bold_terms <- c(
  "negative regulation of alpha-beta T cell differentiation",
  "regulation of alpha-beta T cell differentiation",
  "negative regulation of lymphocyte differentiation",
  "negative regulation of T cell differentiation",
  "myeloid leukocyte differentiation",
  "regulation of lymphocyte differentiation",
  "regulation of leukocyte differentiation",
  "negative regulation of alpha-beta T cell activation",
  "negative regulation of CD4-positive, alpha-beta T cell differentiation",
  "negative regulation of leukocyte differentiation",
  "regulation of alpha-beta T cell activation",
  "negative regulation of lymphocyte activation")

# immune term retained in the motif set (highlighted grey rather than yellow)
grey_term <- "myeloid leukocyte differentiation"

y_faces <- ifelse(
  levels(plot_data$Description) %in% c(bold_terms, grey_term),
  "bold", "plain")

rect_df <- data.frame(
  ypos = which(levels(plot_data$Description) %in% setdiff(bold_terms, grey_term)))

rect_grey_df <- data.frame(
  ypos = which(levels(plot_data$Description) %in% grey_term))

# plot the dotplot
ggplot(plot_data, aes(x = x, y = Description,
                      size = Count, colour = p.adjust)) +
  geom_rect(data = rect_df, inherit.aes = FALSE,
            aes(ymin = ypos - 0.5, ymax = ypos + 0.5,
                xmin = -Inf, xmax = Inf),
            fill = "#FFFF00", alpha = 0.35) +
  geom_rect(data = rect_grey_df, inherit.aes = FALSE,
            aes(ymin = ypos - 0.5, ymax = ypos + 0.5,
                xmin = -Inf, xmax = Inf),
            fill = "#BBBBBB", alpha = 0.35) +
  geom_point() +
  facet_wrap(~ subset, ncol = 2, strip.position = "bottom") +
  scale_colour_gradient(low = "#1a80bb", high = "#cfe3f0",
                        name = "adj. p") +
  scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
  scale_x_continuous(limits = c(0.5, 1.5), breaks = NULL) +
  xlab(NULL) +
  ylab("Top GO Biological Process Terms") +
  ggtitle("HNF1A GO Enrichment (Pancreatic Islets): All Peaks vs Motif-Containing Peaks") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5, margin = margin(b = 25)),
    panel.grid.minor = element_blank(),
    strip.text   = element_text(face = "bold", size = 16),                                # x-axis titles
    axis.title.y = element_text(face = "bold", size = 16),                                # y-axis title
    axis.text.y = element_text(face = y_faces, colour = "black", size = 16),              # GO labels
    legend.title = element_text(face = "bold", size = 14),                                # key titles
    legend.text  = element_text(face = "bold", size = 13))                                # key labels

