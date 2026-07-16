# side-by-side binary (motif vs no motif) boxplots with significance bracket

# load libraries
library(readr)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(patchwork)

# load data
islets <- read_csv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_FIMO_AllPeaks.csv", show_col_types = FALSE)
hepg2  <- read_csv("HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_FIMO_AllPeaks.csv", show_col_types = FALSE)

# rebuild log_peak_score (not stored in the 20-column CSV)
islets <- islets %>% mutate(log_peak_score = log10(peak_score + 1))
hepg2  <- hepg2  %>% mutate(log_peak_score = log10(peak_score + 1))

# correct orientation
islets$motif_group <- factor(islets$motif_group, levels = c("No Motif", "Motif"))
hepg2$motif_group  <- factor(hepg2$motif_group,  levels = c("No Motif", "Motif"))

# reusable builder for one tissue's binary boxplot with sig bracket
build_binary_plot <- function(df, title, box_col, dot_col) {
  
  mw_test <- wilcox.test(peak_score ~ motif_group, data = df, exact = FALSE)
  mw_p    <- signif(mw_test$p.value, 3)
  
  label <- case_when(
    mw_test$p.value < 0.001 ~ "***",
    mw_test$p.value < 0.01  ~ "**",
    mw_test$p.value < 0.05  ~ "*",
    TRUE                    ~ "ns")
  
  y_top  <- max(df$log_peak_score, na.rm = TRUE)
  y_step <- 0.06 * y_top
  y_pos  <- y_top + y_step * (1 + -0.3)
  
  set.seed(42)
  ggplot(df, aes(x = motif_group, y = log_peak_score)) +
    geom_boxplot(outlier.alpha = 0.3, fill = box_col, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1, colour = dot_col) +
    geom_signif(xmin = 1, xmax = 2, y_position = y_pos, annotations = label,
                tip_length = 0.01, textsize = 6, vjust = 0.5, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.04))) +
    theme_classic(base_size = 13) +
    labs(x = "HNF1A Motif Presence",
         y = bquote(bold(log[10]*"(Peak Score + 1)")),
         title = title,
         subtitle = paste0("Mann-Whitney U, p = ", mw_p)) +
    theme(
      plot.title    = element_text(face = "bold", margin = margin(b = 8), size = 20),
      plot.subtitle = element_text(margin = margin(b = 20), size = 19),
      axis.title.x  = element_text(face = "bold", size = 18),
      axis.title.y  = element_text(face = "bold", size = 18),
      axis.text     = element_text(face = "bold", size = 17),
      axis.text.x   = element_text(face = "bold", size = 17))}

# build each panel
plot_islets_binary <- build_binary_plot(
  islets, "Islets: Peak Score By Motif Presence",
  box_col = "#6a98b3", dot_col = "#2c2f33")

plot_hepg2_binary <- build_binary_plot(
  hepg2, "HepG2: Peak Score By Motif Presence",
  box_col = "#c0544a", dot_col = "#2c2f33")

# combine side by side
plot_islets_binary + plot_hepg2_binary +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"))

