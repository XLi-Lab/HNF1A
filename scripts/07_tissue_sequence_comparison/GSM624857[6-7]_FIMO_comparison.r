# this script is to compare the actual HNF1A motif sequences 
# and to see if the matched sequence is different between tissues

# load libraries
library(readr)
library(VennDiagram)
library(grid)
library(dplyr)
library(ggplot2)

# load data
GSM6248576_FIMO <- read_tsv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_MA0046.1.tsv", comment = "#", show_col_types = FALSE)
GSM6248577_FIMO <- read_tsv("HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_MA0046.1.tsv", comment = "#", show_col_types = FALSE)
GSM6248576 <- read_csv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_peakcentric.csv", show_col_types = FALSE)
GSM6248577 <- read_csv("HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_peakcentric.csv", show_col_types = FALSE)

# column name consistency
GSM6248576_FIMO <- GSM6248576_FIMO %>% rename(q_value = `q-value`)
GSM6248577_FIMO <- GSM6248577_FIMO %>% rename(q_value = `q-value`)

# unique matched motif sequences per tissue (FIMO orientation)
islet_seqs <- unique(GSM6248576_FIMO$matched_sequence)
hep_seqs   <- unique(GSM6248577_FIMO$matched_sequence)

# computing the different lengths
length(islet_seqs)
length(hep_seqs)

length(intersect(islet_seqs, hep_seqs))
length(setdiff(islet_seqs, hep_seqs))
length(setdiff(hep_seqs, islet_seqs))

# plotting venn diagram to visualise sequence overlap (shared vs tissue-unique motif sequences)
venn.plot <- VennDiagram::venn.diagram(
  x = list(
    Islet = islet_seqs,
    HepG2 = hep_seqs),
  category.names = c("Islet", "HepG2"),
  filename = NULL,
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  margin = 0.1)
grid::grid.newpage()
grid::grid.draw(venn.plot)

# building the peak --> gene lookup
# first make peak IDs match FIMO sequence_name format
islet_peak_annot <- GSM6248576 %>%
  mutate(sequence_name = sub(":", "_", peak_id)) %>%
  select(sequence_name, nearest_proteincoding_gene)

hep_peak_annot <- GSM6248577 %>%
  mutate(sequence_name = sub(":", "_", peak_id)) %>%
  select(sequence_name, nearest_proteincoding_gene)

# then join every FIMO motif hit to nearest protein-coding genes
islet_fimo_gene <- GSM6248576_FIMO %>%
  left_join(islet_peak_annot, by = "sequence_name") %>%
  mutate(tissue = "Islet")

hep_fimo_gene <- GSM6248577_FIMO %>%
  left_join(hep_peak_annot, by = "sequence_name") %>%
  mutate(tissue = "HepG2")

# define gene groups using ONLY motif-positive peaks
islet_genes <- islet_fimo_gene %>%
  pull(nearest_proteincoding_gene) %>%
  na.omit() %>%
  unique()

hep_genes <- hep_fimo_gene %>%
  pull(nearest_proteincoding_gene) %>%
  na.omit() %>%
  unique()

# computing the different lengths
shared_genes <- intersect(islet_genes, hep_genes)
islet_only   <- setdiff(islet_genes, hep_genes)
hep_only     <- setdiff(hep_genes, islet_genes)

length(islet_genes)
length(hep_genes)
length(shared_genes)
length(islet_only)
length(hep_only)

# consensus sequence
consensus <- "GGTTAATNATTAAC"

# helper function computes hamming distance of a sequence to the consensus
hamming_to_consensus <- function(seqs, consensus) {
  cons <- strsplit(consensus, "")[[1]]
  vapply(seqs, function(s) {
    s <- strsplit(toupper(s), "")[[1]]
    if (length(s) != length(cons)) return(NA_real_)
    sum(cons != "N" & s != cons)}, numeric(1))}

# create combined FIMO dataframe
GSM6248576.7_Combined_FIMO <- bind_rows(islet_fimo_gene, hep_fimo_gene) %>%
  mutate(
    target_group = case_when(
      nearest_proteincoding_gene %in% shared_genes ~ "Shared",
      tissue == "Islet" & nearest_proteincoding_gene %in% islet_only  ~ "Islet-Specific",
      tissue == "HepG2" & nearest_proteincoding_gene %in% hep_only    ~ "HepG2-Specific",
      TRUE ~ NA_character_),
    target_group = factor(target_group,
                          levels = c("Islet-Specific", "Shared", "HepG2-Specific")),
    hamming_distance = hamming_to_consensus(matched_sequence, consensus),
    score   = as.numeric(score),
    q_value = as.numeric(q_value)) %>%
  filter(!is.na(target_group), !is.na(hamming_distance))

# check group size
table(GSM6248576.7_Combined_FIMO$tissue, GSM6248576.7_Combined_FIMO$target_group)
table(GSM6248576.7_Combined_FIMO$target_group)

# plot hamming distance for each peak in islet-specific, shared and hepatocyte-specific
set.seed(42)
ggplot(GSM6248576.7_Combined_FIMO, aes(x = target_group, y = hamming_distance)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.25, size = 1) +
  theme_classic() +
  labs(
    x = "Target Gene Group",
    y = "Hamming Distance From HNF1A Consensus",
    title = "HNF1A Motif Divergence In Tissue-Specific & Shared Targets")

# statistical test to see if hamming distance difference between the 3 groups is significant or not
kruskal.test(hamming_distance ~ target_group, data = GSM6248576.7_Combined_FIMO)

pairwise.wilcox.test(
  GSM6248576.7_Combined_FIMO$hamming_distance,
  GSM6248576.7_Combined_FIMO$target_group,
  p.adjust.method = "BH")

# save output
write.csv(GSM6248576.7_Combined_FIMO, "GSM624857[6-7]_Outputs/GSM624857[6-7]_Tissue_Sequence_Comparison.csv", row.names = FALSE)
