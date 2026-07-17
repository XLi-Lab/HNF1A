# this r script computes the genomic-feature distribution of HepG2 HNF1A peaks (GSM6248577),
# tests enrichment against a genome background: observed vs expected per category
# output saved to GSM6248577_Enrichment_Profile.csv

# load libraries
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(dplyr)
library(tidyr)

# count peaks per genomic category (observed distribution)
cats <- c("Promoter","5UTR","3UTR","Exon","Intron","Downstream", "Upstream", "Intergenic")

# slightly different code to GSM6248576 due to the presence of NAs
obs <- table(GSM6248577$peak_overlap_custom, useNA = "ifany")
obs <- obs[cats]            # enforces ordering
obs[is.na(obs)] <- 0        # if a category is absent in this dataset
stopifnot(sum(obs) + sum(is.na(GSM6248577$peak_overlap_custom)) == nrow(GSM6248577))

reduce <- GenomicRanges::reduce

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# use standard autosomal + sex chromosomes only
keep <- grep("^chr([0-9]+|X|Y)$", seqlevels(txdb), value = TRUE)
txdb2 <- keepSeqlevels(txdb, keep, pruning.mode = "coarse")

# genome size (sum of seqlengths across standard chromosomes)
genome_size <- sum(seqlengths(txdb2), na.rm = TRUE)

# gene models (strand-aware)
genes0 <- genes(txdb2)

# promoters which are already defined (around TSS)
prom_gr <- reduce(promoters(genes0, upstream = 2000, downstream = 1000))

# upstream region found -3000 to -2001 relative to TSS
upstream_gr <- reduce(setdiff(promoters(genes0, upstream = 3000, downstream = 0), promoters(genes0, upstream = 2000, downstream = 0)))

# downstream region found 0..+3000 from TES
down_gr <- reduce(flank(genes0, width = 3000, start = FALSE))

# exons
exon_gr <- reduce(exons(txdb2))

# UTRs
utr5_gr <- tryCatch(reduce(unlist(fiveUTRsByTranscript(txdb2), use.names = FALSE)),
                    error = function(e) GRanges())
utr3_gr <- tryCatch(reduce(unlist(threeUTRsByTranscript(txdb2), use.names = FALSE)),
                    error = function(e) GRanges())

# genes (full gene body)
gene_gr <- reduce(genes(txdb2))

# introns = gene body minus exons
intron_gr <- reduce(setdiff(gene_gr, exon_gr))

# make disjoint bins to avoid double counting, with a priority scheme
# priority: promoter > UTRs > exon > intron > downstream > intergenic
utr5_only   <- reduce(setdiff(utr5_gr, prom_gr))
utr3_only   <- reduce(setdiff(utr3_gr, prom_gr))
exon_only   <- reduce(setdiff(exon_gr, union(prom_gr, union(utr5_only, utr3_only))))
intron_only <- reduce(setdiff(intron_gr, union(prom_gr, union(utr5_only, union(utr3_only, exon_only)))))
upstream_only <- reduce(setdiff(upstream_gr, union(prom_gr, union(utr5_only, union(utr3_only, union(exon_only, intron_only))))))
down_only <- reduce(setdiff(down_gr, union(prom_gr, union(utr5_only, union(utr3_only, union(exon_only, union(intron_only, upstream_only)))))))

# everything else becomes intergenic
covered <- reduce(union(prom_gr, union(utr5_only, union(utr3_only, union(exon_only, union(intron_only, union(upstream_only, down_only)))))))
intergenic_bp <- genome_size - sum(width(covered))

exp_bp <- c(
  Promoter   = sum(width(prom_gr)),
  `5UTR`     = sum(width(utr5_only)),
  `3UTR`     = sum(width(utr3_only)),
  Exon       = sum(width(exon_only)),
  Intron     = sum(width(intron_only)),
  Downstream = sum(width(down_only)),
  Upstream   = sum(width(upstream_only)),
  Intergenic = intergenic_bp)

# align expected to observed to prevent NA
exp_bp_aligned <- exp_bp[cats]
stopifnot(!any(is.na(exp_bp_aligned)))

# converts basepair counts into genome percentage
# what fraction of the genome (in bp) falls into each category?
round(100 * exp_bp_aligned / sum(exp_bp_aligned), 3)

## global chi-square test
exp_prob <- exp_bp_aligned / sum(exp_bp_aligned)
chisq_res <- chisq.test(x = as.numeric(obs), p = exp_prob)

## odds ratio per category
obs_total <- sum(obs)

exp_prob <- exp_bp_aligned / sum(exp_bp_aligned)
exp_counts <- round(obs_total * exp_prob)

# fix rounding so totals match exactly
diff_ct <- obs_total - sum(exp_counts)
if (diff_ct != 0) {
  exp_counts[which.max(exp_counts)] <- exp_counts[which.max(exp_counts)] + diff_ct}

# quick checks
exp_counts
sum(exp_counts)       # should equal obs_total

# computing per-category enrichment statistics
# asking is this category enriched or depleted in peaks compared to what we'd expect from the genome background?

# create a contingency table to compare observed vs expected
calc_or <- function(cat) {
  m <- matrix(c(
    obs[cat],          obs_total - obs[cat],
    exp_counts[cat],   obs_total - exp_counts[cat]
  ), nrow = 2, byrow = TRUE)
  
  # run Fisher's exact test
  ft <- fisher.test(m)
  
  data.frame(
    category = cat,
    OR = unname(ft$estimate),
    p_value = ft$p.value,
    CI_low = ft$conf.int[1],
    CI_high = ft$conf.int[2])}

category_enrichment_results <- do.call(rbind, lapply(names(obs), calc_or))
category_enrichment_results$padj <- p.adjust(category_enrichment_results$p_value, method = "BH")
category_enrichment_results[order(category_enrichment_results$padj), ]

# clean summary table of your enrichment analysis
data.frame(
  category = names(obs),
  obs_n = as.integer(obs),
  obs_pct = round(100 * obs / sum(obs), 2),
  exp_pct = round(100 * exp_prob, 2),
  exp_n = as.integer(exp_counts))

# plotting
category_distribution_df <- data.frame(
  category = names(obs),
  Observed = as.numeric(obs) / sum(obs) * 100,
  Expected = exp_prob * 100)

category_distribution_long_df <- tidyr::pivot_longer(
  category_distribution_df,
  cols = c("Observed", "Expected"),
  names_to = "Type",
  values_to = "Percent")

ggplot(category_distribution_long_df, aes(x = category, y = Percent, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Percentage") +
  xlab("Genomic Category") +
  theme_minimal(base_size = 14)

# fold enrichment
fold_enrichment <- (as.numeric(obs) / sum(obs)) / exp_prob

fe_df <- data.frame(
  category = names(obs),
  obs_n    = as.integer(obs),
  FoldEnrichment = (as.numeric(obs) / sum(obs)) / exp_prob) |>
  mutate(category = reorder(category, FoldEnrichment))

ggplot(fe_df, aes(x = category, y = FoldEnrichment)) +
  geom_col(width = 0.7, fill = "#D85A30") +               
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = paste0("n=", obs_n)),
            hjust = -0.15, size = 3.1, colour = "grey30") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "HepG2 Peaks",
       x = NULL, y = "Fold Enrichment (Observed / Expected)") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

# odds ratio forest
ggplot(category_enrichment_results, aes(x = category, y = OR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  ylab("Odds Ratio") +
  theme_minimal(base_size = 14)

# build a distance to TSS plot
GSM6248577$peak_overlap_custom <- factor(
  GSM6248577$peak_overlap_custom,
  levels = c("Promoter","5UTR","3UTR","Exon","Intron","Downstream","Upstream","Intergenic"))

ggplot(
  subset(GSM6248577,
         !is.na(peak_distanceToTSS) &
           peak_distanceToTSS >= -10000 &
           peak_distanceToTSS <= 10000),
  aes(x = peak_distanceToTSS,
      y = peak_overlap_custom,
      colour = peak_overlap_custom)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Distance to TSS (bp)",
       y = "Annotation",
       colour = "Annotation") +
  theme_minimal(base_size = 14)

# save output
enrichment_out <- data.frame(
  category = names(obs),
  obs      = as.integer(obs),
  exp_prob = as.numeric(exp_prob))

write.csv(enrichment_out, "HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_Enrichment_Profile.csv", row.names = FALSE)
