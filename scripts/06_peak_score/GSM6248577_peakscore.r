# this r script tests HepG2 HNF1A ChIP-Seq peak score against motif presence and motif count

# load libraries
library(dplyr)
library(ggplot2)
library(readr)
library(FSA)

# read the per-peak motif table written by GSM6248577_FIMO.r
all_peaks_with_counts <- read_csv("HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_FIMO_AllPeaks.csv", show_col_types = FALSE)

# restore factor ordering (as CSV stores motif_group as plain text)
all_peaks_with_counts <- all_peaks_with_counts %>%
  mutate(motif_group = factor(motif_group, levels = c("No Motif", "Motif")))

# log transform peak score for plotting (skewed data)
all_peaks_with_counts <- all_peaks_with_counts %>%
  mutate(log_peak_score = log10(peak_score + 1))

# test1, BINARY (motif vs no motif)
# Mann-Whitney U
mw_test <- wilcox.test(peak_score ~ motif_group,
                       data = all_peaks_with_counts,
                       exact = FALSE)
print(mw_test)
mw_p <- signif(mw_test$p.value, 3)

set.seed(42)
plot_hepg2_binary <- ggplot(all_peaks_with_counts, aes(x = motif_group, y = log_peak_score)) +
  geom_boxplot(outlier.alpha = 0.3, fill = "#c9cdd6", width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1, colour = "#2c2f33") +
  theme_classic(base_size = 13) +
  labs(x = "HNF1A Motif Presence",
       y = bquote(bold(log[10]*"(Peak Score + 1)")),
       title = "GSM6248577 (HepG2): Peak Score By Motif Presence",
       subtitle = paste0("Mann-Whitney U, p = ", mw_p)) +
  theme(
    plot.title    = element_text(face = "bold", margin = margin(b = 8)),
    plot.subtitle = element_text(margin = margin(b = 20)),
    axis.title.x  = element_text(face = "bold", size = 14),
    axis.title.y  = element_text(face = "bold", size = 14))

plot_hepg2_binary

# test2, COUNT (number of motifs per peak)
# bin motif counts (HepG2 uses finer bins than islets: 0/1/2/3/4+)
all_peaks_with_counts <- all_peaks_with_counts %>%
  mutate(motif_count_group = case_when(
    distinct_motifs == 0 ~ "0",
    distinct_motifs == 1 ~ "1",
    distinct_motifs == 2 ~ "2",
    distinct_motifs == 3 ~ "3",
    distinct_motifs >= 4 ~ "4+"),
    motif_count_group = factor(motif_count_group,
                               levels = c("0", "1", "2", "3", "4+")))

# check bin sizes first
table(all_peaks_with_counts$motif_count_group)

# Kruskal-Wallis + Dunn's
kw_test <- kruskal.test(peak_score ~ motif_count_group, data = all_peaks_with_counts)
print(kw_test)

# as KW test is significant, do Dunn's test
dunn_test <- dunnTest(peak_score ~ motif_count_group,
                      data = all_peaks_with_counts, method = "bh")
print(dunn_test)

# Spearman correlation (compute before plotting so subtitle works)
sp_test <- cor.test(all_peaks_with_counts$distinct_motifs,
                    all_peaks_with_counts$peak_score,
                    method = "spearman", exact = FALSE)
print(sp_test)

# annotated grouped boxplot
kw_p   <- signif(kw_test$p.value, 3)
sp_rho <- signif(sp_test$estimate, 3)
sp_p   <- signif(sp_test$p.value, 3)

set.seed(42)
ggplot(all_peaks_with_counts, aes(x = motif_count_group, y = log_peak_score)) +
  geom_boxplot(outlier.alpha = 0.3, fill = "#c9cdd6") +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1, colour = "#2c2f33") +
  theme_classic(base_size = 13) +
  labs(x = "Number of HNF1A Motifs Per Peak",
       y = bquote(bold(log[10]*"(Peak Score + 1)")),
       title = "GSM6248577 (HepG2): Peak Score By Motif Count",
       subtitle = paste0("Kruskal-Wallis p = ", kw_p,
                         "  |  Spearman ρ = ", sp_rho, ", p = ", sp_p)) +
  theme(
    plot.title    = element_text(face = "bold", margin = margin(b = 8), size = 20),
    plot.subtitle = element_text(margin = margin(b = 20), size = 16),
    axis.title.x  = element_text(face = "bold", size = 16),
    axis.title.y  = element_text(face = "bold", size = 16),
    axis.text     = element_text(face = "bold", size = 13))

# extra: group sizes and medians
all_peaks_with_counts %>%
  group_by(motif_group) %>%
  summarise(n = n(),
            median_score = median(peak_score),
            median_log_score = median(log_peak_score))

