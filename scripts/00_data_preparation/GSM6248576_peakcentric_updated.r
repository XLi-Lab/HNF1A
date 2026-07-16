# the r script builds the master one-row-per-peak table for islet HNF1A peaks (GSM6248576)
# it reads the raw peak BED, annotates each peak's genomic location,
# assigns nearest gene and nearest protein-coding gene with TSS distances,
# lists genes within ±10 kb, and carries the peak score
# output is saved as GSM6248576_peakcentric.csv

# load libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)  
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(EnsDb.Hsapiens.v86)})

# load TxDb database
# load ENSEMBL database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
edb  <- EnsDb.Hsapiens.v86

# load peaks
# reads BED file into a GRanges
# keeps only standard chromosomes
# sorts chromosomes first and then by peak ranges
peaks <- readPeakFile("/Users/Manu/Downloads/GSE206240_RAW/HNF1A/GSM6248576_Islets_HNF1A_ab96777_peaks.bed")
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- dropSeqlevels(peaks, "chrM", pruning.mode = "coarse")
peaks <- sortSeqlevels(peaks)
peaks <- sort(peaks)

# create a 10,000 base pair window around each peak midpoint
WINDOW <- 10000

# peak_id
peak_id <- paste0(as.character(seqnames(peaks)), ":", start(peaks), "-", end(peaks))

# midpoint (1bp)
mid <- start(peaks) + floor((width(peaks) - 1) / 2)
peak_mid <- GRanges(
  seqnames = seqnames(peaks),
  ranges   = IRanges(start = mid, width = 1),
  peak_id  = peak_id)

# +/-10kb window around midpoint
peak_win10 <- resize(peak_mid, width = 2 * WINDOW + 1, fix = "center")

# find nereast protein-coding gene using Ensembl database
# TxDb.Hsapiens.UCSC.hg38.knownGene object has coordinates, but it does not tell you gene biotype, Ensembl does.
proteincoding_genes <- genes(
  edb,
  filter = AnnotationFilter::GeneBiotypeFilter("protein_coding"),
  return.type = "GRanges")

# make the chromosome style consistent between peak_mid (chr1, chr2, chr3...) and proteincoding_genes (1, 2, 3...)
GenomeInfoDb::seqlevelsStyle(proteincoding_genes) <- GenomeInfoDb::seqlevelsStyle(peaks)

# keep only standard chromosomes to match peaks and sort
proteincoding_genes <- keepStandardChromosomes(proteincoding_genes, pruning.mode = "coarse")
proteincoding_genes <- sortSeqlevels(proteincoding_genes)
proteincoding_genes <- sort(proteincoding_genes)

# map Ensembl gene IDs to gene symbols
proteincoding_symbols <- mapIds(
  org.Hs.eg.db,
  keys = unique(mcols(proteincoding_genes)$gene_id),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first")

# find nearest protein coding gene to each peak midpoint (based on gene body)
nearestproteincoding_gene <- nearest(peak_mid, proteincoding_genes, ignore.strand = TRUE)

# extract Ensembl ID and symbol of that nearest gene
nearestproteincoding_ensembl <- mcols(proteincoding_genes)$gene_id[nearestproteincoding_gene]
nearestproteincoding_symbol  <- unname(proteincoding_symbols[nearestproteincoding_ensembl])

# build TSS positions for protein-coding genes, get the TSS for the nearest protein-coding gene, compute midpoint-to-TSS distance
proteincoding_tss <- promoters(proteincoding_genes, upstream = 0, downstream = 1)
nearestproteincoding_tss <- proteincoding_tss[nearestproteincoding_gene]
nearestproteincoding_distance <- start(peak_mid) - start(nearestproteincoding_tss)

# fix the sign for minus-strand genes
minus <- as.character(strand(nearestproteincoding_tss)) == "-"
nearestproteincoding_distance[minus] <- start(nearestproteincoding_tss)[minus] - start(peak_mid)[minus]

# get gene models
genes_gr <- genes(txdb)
genes_gr$ENTREZID <- as.character(genes_gr$gene_id)

# TSS points (1bp) for each gene
tss_gr <- promoters(genes_gr, upstream = 0, downstream = 1)

# Entrez -> SYMBOL
entrez_to_symbol <- function(entrez_vec) {
  entrez_vec <- unique(as.character(entrez_vec))
  entrez_vec <- entrez_vec[!is.na(entrez_vec) & entrez_vec != ""]
  if (length(entrez_vec) == 0) return(character(0))
  syms <- mapIds(
    org.Hs.eg.db,
    keys = entrez_vec,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first")
  unique(na.omit(unname(syms)))}

# a) genes whose TSS is within +/-10kb of peak midpoint
hits_tss <- findOverlaps(peak_win10, tss_gr, ignore.strand = TRUE)
tss_list_entrez <- split(tss_gr$ENTREZID[subjectHits(hits_tss)], queryHits(hits_tss))

# b) genes whose gene body overlaps the +/-10kb window
hits_body <- findOverlaps(peak_win10, genes_gr, ignore.strand = TRUE)
body_list_entrez <- split(genes_gr$ENTREZID[subjectHits(hits_body)], queryHits(hits_body))

# c) nearest gene to each peak midpoint
nearest_idx <- nearest(peak_mid, genes_gr, ignore.strand = TRUE)
nearest_entrez <- genes_gr$ENTREZID[nearest_idx]
nearest_sym_map <- mapIds(
  org.Hs.eg.db,
  keys = unique(nearest_entrez),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first")

# build one-row-per-peak output (main table)
GSM6248576 <- data.frame(
  peak_id = peak_id,
  peak_chr = as.character(seqnames(peaks)),
  peak_start = start(peaks),
  peak_end = end(peaks),
  peak_mid = start(peak_mid),
  peak_score = mcols(peaks)$V5,
  win10_start = start(peak_win10),
  win10_end = end(peak_win10),
  genes_TSS_10kb = "",
  genes_geneBodyOverlap_10kb = "",
  nearest_gene = as.character(unname(nearest_sym_map[nearest_entrez])),
  nearest_proteincoding_gene = nearestproteincoding_symbol,
  nearest_proteincoding_gene_distance = nearestproteincoding_distance,
  stringsAsFactors = FALSE)

# fill gene lists per peak index (queryHits are indices into peak_win10)
for (i in seq_len(nrow(GSM6248576))) {
  if (!is.null(tss_list_entrez[[as.character(i)]])) {
    GSM6248576$genes_TSS_10kb[i] <- paste(entrez_to_symbol(tss_list_entrez[[as.character(i)]]), collapse = ",")}
  if (!is.null(body_list_entrez[[as.character(i)]])) {
    GSM6248576$genes_geneBodyOverlap_10kb[i] <- paste(entrez_to_symbol(body_list_entrez[[as.character(i)]]), collapse = ",")}}

# peak overlap annotation: where the PEAK lies in genome based on TxDb hg38
peakAnno <- annotatePeak(
  peaks,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db",
  tssRegion = c(-2000, 1000))  # promoter definition

pa <- as.data.frame(peakAnno)
pa$peak_id <- paste0(pa$seqnames, ":", pa$start, "-", pa$end)

# join onto GSM6248576 by peak_id
GSM6248576 <- merge(
  GSM6248576,
  pa[, c("peak_id", "annotation", "distanceToTSS")],
  by = "peak_id",
  all.x = TRUE,
  sort = FALSE)

names(GSM6248576)[names(GSM6248576) == "annotation"] <- "peak_overlap_type"
names(GSM6248576)[names(GSM6248576) == "distanceToTSS"] <- "peak_distanceToTSS"

# simplified overlap category for summary statistics
GSM6248576$peak_overlap_simple <- dplyr::case_when(
  is.na(GSM6248576$peak_overlap_type)                  ~ NA_character_,
  grepl("Promoter",  GSM6248576$peak_overlap_type)    ~ "Promoter",
  grepl("5' UTR",    GSM6248576$peak_overlap_type)    ~ "5UTR",
  grepl("3' UTR",    GSM6248576$peak_overlap_type)    ~ "3UTR",
  grepl("Exon",      GSM6248576$peak_overlap_type)    ~ "Exon",
  grepl("Intron",    GSM6248576$peak_overlap_type)    ~ "Intron",
  grepl("Downstream",GSM6248576$peak_overlap_type)    ~ "Downstream",
  grepl("Intergenic",GSM6248576$peak_overlap_type)    ~ "Intergenic",
  TRUE                                                 ~ "Other")

# match peak names for FIMO analysis
GSM6248576 <- GSM6248576 %>%
  mutate(sequence_name = gsub(":", "_", peak_id))

# add upstream bin (subset of intergenic peaks upstream of promoter)
GSM6248576$peak_overlap_custom <- dplyr::case_when(
  is.na(GSM6248576$peak_overlap_simple) ~ NA_character_,
  
  # keep existing annotations
  GSM6248576$peak_overlap_simple != "Intergenic" ~ GSM6248576$peak_overlap_simple,
  
  # redefine subset of intergenic peaks as upstream of promoter
  GSM6248576$peak_overlap_simple == "Intergenic" &
    !is.na(GSM6248576$peak_distanceToTSS) &
    GSM6248576$peak_distanceToTSS >= -3000 &
    GSM6248576$peak_distanceToTSS < -2000 ~ "Upstream",
  
  # remaining intergenic peaks stay intergenic
  TRUE ~ "Intergenic")

# sort via order that makes sense
GSM6248576$peak_overlap_custom <- factor(
  GSM6248576$peak_overlap_custom,
  levels = c(
    "Promoter",
    "5UTR",
    "3UTR",
    "Exon",
    "Intron",
    "Upstream",
    "Downstream",
    "Intergenic"))

# save output
write.csv(GSM6248576, "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_peakcentric.csv", row.names = FALSE)

