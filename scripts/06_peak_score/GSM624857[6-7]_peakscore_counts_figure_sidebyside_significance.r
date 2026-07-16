# side-by-side boxplots with significance bracket for motif-containg and motif > 1
# plus Spearman dose-response comparison (all peaks vs motif-only), peak score vs motif count

# load libraries
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(FSA)

# load data
islets <- read_csv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_FIMO_AllPeaks.csv", show_col_types = FALSE)
hepg2  <- read_csv("HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_FIMO_AllPeaks.csv", show_col_types = FALSE)

# rebuild columns peakscore derives (not stored in the 20-column CSV)
islets <- islets %>%
  mutate(log_peak_score = log10(peak_score + 1),
         motif_count_group = case_when(
           distinct_motifs == 0 ~ "0",
           distinct_motifs == 1 ~ "1",
           distinct_motifs == 2 ~ "2",
           distinct_motifs >= 3 ~ "3+"))

hepg2 <- hepg2 %>%
  mutate(log_peak_score = log10(peak_score + 1),
         motif_count_group = case_when(
           distinct_motifs == 0 ~ "0",
           distinct_motifs == 1 ~ "1",
           distinct_motifs == 2 ~ "2",
           distinct_motifs == 3 ~ "3",
           distinct_motifs >= 4 ~ "4+"))

# enforce correct factor level order per tissue
islets$motif_count_group <- factor(islets$motif_count_group, levels = c("0", "1", "2", "3+"))
hepg2$motif_count_group  <- factor(hepg2$motif_count_group,  levels = c("0", "1", "2", "3", "4+"))

# reusable builder for one tissue's grouped boxplot with sig brackets
build_grouped_plot <- function(df, title, box_col, dot_col, only_vs_zero) {
  
  kw_test   <- kruskal.test(peak_score ~ motif_count_group, data = df)
  dunn_test <- dunnTest(peak_score ~ motif_count_group, data = df, method = "bh")
  
  lv   <- levels(df$motif_count_group)
  xpos <- setNames(seq_along(lv), lv)
  
  sig_brackets <- dunn_test$res %>%
    mutate(Comparison = str_replace_all(Comparison, " ", "")) %>%
    separate(Comparison, into = c("g1", "g2"), sep = "-", remove = FALSE) %>%
    filter(P.adj < 0.05)
  
  if (only_vs_zero) sig_brackets <- filter(sig_brackets, g1 == "0" | g2 == "0")
  
  sig_brackets <- sig_brackets %>%
    mutate(
      xmin  = as.numeric(xpos[g1]),
      xmax  = as.numeric(xpos[g2]),
      span  = abs(xmax - xmin),
      label = case_when(
        P.adj < 0.001 ~ "***",
        P.adj < 0.01  ~ "**",
        P.adj < 0.05  ~ "*",
        TRUE          ~ "ns")) %>%
    arrange(span, xmin)
  
  y_top  <- max(df$log_peak_score, na.rm = TRUE)
  y_step <- 0.06 * y_top
  sig_brackets$y <- y_top + y_step * (seq_len(nrow(sig_brackets)) + -0.3)
  
  kw_p <- signif(kw_test$p.value, 3)
  
  set.seed(42)
  p <- ggplot(df, aes(x = motif_count_group, y = log_peak_score)) +
    geom_boxplot(outlier.alpha = 0.3, fill = box_col) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1, colour = dot_col) +
    geom_signif(xmin = sig_brackets$xmin, xmax = sig_brackets$xmax,
                y_position = sig_brackets$y, annotations = sig_brackets$label,
                tip_length = 0.01, textsize = 6, vjust = 0.5, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.04))) +
    theme_classic(base_size = 13) +
    labs(x = "Number of HNF1A Motifs Per Peak",
         y = bquote(bold(log[10]*"(Peak Score + 1)")),
         title = title,
         subtitle = paste0("Kruskal-Wallis p = ", kw_p)) +
    theme(
      plot.title    = element_text(face = "bold", margin = margin(b = 8), size = 20),
      plot.subtitle = element_text(margin = margin(b = 20), size = 19),
      axis.title.x  = element_text(face = "bold", size = 18),
      axis.title.y  = element_text(face = "bold", size = 18),
      axis.text     = element_text(face = "bold", size = 17),
      axis.text.x   = element_text(face = "bold", size = 17))
  
  list(plot = p, dunn = dunn_test$res)}

# build each panel
islets_out <- build_grouped_plot(
  islets, "Islets: Peak Score By Motif Count",
  box_col = "#6a98b3", dot_col = "#2c2f33", only_vs_zero = TRUE)

hepg2_out <- build_grouped_plot(
  hepg2, "HepG2: Peak Score By Motif Count",
  box_col = "#c0544a", dot_col = "#2c2f33", only_vs_zero = FALSE)

plot_islets_grouped <- islets_out$plot
plot_hepg2_grouped  <- hepg2_out$plot

# assemble panels c) and d)
plot_islets_grouped + plot_hepg2_grouped +
  plot_annotation(tag_levels = list(c("c", "d")), tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold"))

# stats, doing Spearman on the full set (incl. 0) and on motif-only peaks (>= 1)
spearman_compare <- function(df, label) {
  full <- cor.test(df$distinct_motifs, df$peak_score,
                   method = "spearman", exact = FALSE)
  
  motif_only <- df %>% filter(distinct_motifs >= 1)
  mo <- cor.test(motif_only$distinct_motifs, motif_only$peak_score,
                 method = "spearman", exact = FALSE)
  
  tibble(
    tissue  = label,
    set     = c("all peaks (incl. 0)", "motif-only (>= 1)"),
    n       = c(nrow(df), nrow(motif_only)),
    rho     = signif(c(full$estimate, mo$estimate), 3),
    p_value = signif(c(full$p.value,  mo$p.value),  3))}

comparison <- bind_rows(
  spearman_compare(islets, "Islets (GSM6248576)"),
  spearman_compare(hepg2,  "HepG2 (GSM6248577)"))

print(comparison)

# count spread among motif-containing peaks (context for the >= 1 result)
table(islets$distinct_motifs[islets$distinct_motifs >= 1])
table(hepg2$distinct_motifs[hepg2$distinct_motifs  >= 1])

