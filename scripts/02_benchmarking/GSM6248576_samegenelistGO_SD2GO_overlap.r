# this r script is about doing GO analysis using the gene list the paper itself provided, 
# and seeing if we can recreate the GO analysis the paper published to set a benchmark

# load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(readxl)
library(VennDiagram)
library(grid)

# load GSM6248576_SD2 file into R as a dataframe
GSM6248576_SD2 <- readxl::read_xlsx("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_SD2_GeneList.xlsx")

# load GSM6248576_SD2_GO file into R as a dataframe
GSM6248576_SD2_GO <- readxl::read_xlsx("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_SD2_GO.xlsx")

# get genes (as ENTREZ IDs)
GSM6248576_SD2_GeneList <- unique(as.character(na.omit(GSM6248576_SD2$`Gene ID`)))

# filter to remove dead/retired IDs that won't behave in downstream enrichment
# keeping genes that actually map
GSM6248576_SD2_GeneList_Mapped <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = GSM6248576_SD2_GeneList,
  keytype = "ENTREZID",
  columns = "SYMBOL")

GSM6248576_SD2_GeneList_MappedFinal <- GSM6248576_SD2_GeneList_Mapped$ENTREZID[!is.na(GSM6248576_SD2_GeneList_Mapped$SYMBOL)]

# conduct GO analysis, this is the GO analysis I AM doing using the SAME gene list as the paper
GSM6248576_SD2_GO_Recreation_SameGeneList <- enrichGO(
  gene          = GSM6248576_SD2_GeneList_MappedFinal,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

# comparison between MY GO and PAPER GO
recreation_paper_GO <- as.data.frame(GSM6248576_SD2_GO_Recreation_SameGeneList)
paper_GO <- GSM6248576_SD2_GO

# extract unique GO IDs
recreation_paper_ids <- unique(na.omit(recreation_paper_GO$ID))
paper_ids <- unique(na.omit(paper_GO$ID))

# terms unique to each
recreation_paper_only <- setdiff(recreation_paper_ids, paper_ids)
paper_only <- setdiff(paper_ids, recreation_paper_ids)
shared_ids <- intersect(recreation_paper_ids, paper_ids)

# make the venn diagram
venn.plot <- venn.diagram(
  x = list(
    Recreation_GO = recreation_paper_ids,
    Published_GO = paper_ids),
  category.names = c("", ""),
  main = "GO Term Concordance: Pipeline Recreation vs. Published (Islets)",
  main.cex = 1.3,
  main.fontface = "bold",
  main.pos = c(0.5, 0.9),
  margin = 0.05,
  filename = NULL,
  fill = c("skyblue", "pink"),
  alpha = 0.5,
  cex = 1.5,
  disable.logging = TRUE)

grid.newpage()
grid.draw(venn.plot)
grid.text("Recreation of Published GO (Same Gene List)",
          x = 0.25, y = 0.02, gp = gpar(fontsize = 12))
grid.text("Published GO",
          x = 0.88, y = 0.02, gp = gpar(fontsize = 12))

# how much did the recreation GO analysis recover compared to the paper?
recovery_islets <- length(shared_ids) / length(paper_ids) * 100
recovery_islets

# save output as CSVs to use for (arm 1, SameGeneListGO) figure creation
write.csv(recreation_paper_GO, "HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_SameGeneListGO.csv", row.names = FALSE)

