# this r script builds the arm-2 benchmarking figure (project-derived protein-coding gene list)
# side-by-side GO term concordance Venns for islets and HepG2

# load libraries
library(VennDiagram)
library(grid)
library(ggplotify)
library(patchwork)
library(ggplot2)
library(readxl)
library(readr)

# read my protein-coding GO (from the per-tissue arm-2 scripts) and the published paper GO
islet_myGO  <- read_csv("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_ProteinCodingGO_AllPeaks.csv", show_col_types = FALSE)
islet_paper <- readxl::read_xlsx("HNF1A_Islets_GSM6248576/GSM6248576_Outputs/GSM6248576_SD2_GO.xlsx")

hep_myGO  <- read_csv("HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_ProteinCodingGO_AllPeaks.csv", show_col_types = FALSE)
hep_paper <- readxl::read_xlsx("HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/GSM6248577_SD2_GO.xlsx")

# derive ID sets
islet_myGO_ids  <- unique(na.omit(trimws(islet_myGO$ID)))
islet_paper_ids <- unique(na.omit(trimws(islet_paper$ID)))
hep_myGO_ids    <- unique(na.omit(trimws(hep_myGO$ID)))
hep_paper_ids   <- unique(na.omit(trimws(hep_paper$ID)))

# build venn diagrams (one for Islets and one for HepG2)
venn_islets <- venn.diagram(
  x = list(
    My_GO     = islet_myGO_ids,
    Published = islet_paper_ids),
  category.names  = c("", ""),
  filename        = NULL,
  fill            = c("skyblue", "salmon"),
  alpha           = 0.5,
  cex             = 1.6,
  fontface        = "bold",
  margin          = 0.05,
  disable.logging = TRUE)

venn_hepg2 <- venn.diagram(
  x = list(
    My_GO     = hep_myGO_ids,
    Published = hep_paper_ids),
  category.names  = c("", ""),
  filename        = NULL,
  fill            = c("skyblue", "salmon"),
  alpha           = 0.5,
  cex             = 1.6,
  fontface        = "bold",
  margin          = 0.05,
  disable.logging = TRUE)

# convert to ggplot objects
gg_islets <- as.ggplot(grobTree(venn_islets)) +
  labs(title = "Islets (GSM6248576)") +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15,
                              margin = margin(b = -20)))

gg_hepg2 <- as.ggplot(grobTree(venn_hepg2)) +
  labs(title = "HepG2 (GSM6248577)") +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15,
                              margin = margin(b = -20)))

# add bottom labels (same left/right order for both panels)
gg_islets <- gg_islets +
  annotate("text", x = 0.18, y = 0.1,
           label = "Present Analysis\n(Protein-Coding Genes)",
           fontface = "bold", size = 4.5) +
  annotate("text", x = 0.82, y = 0.1,
           label = "Published\n(Ng et al., 2024)",
           fontface = "bold", size = 4.5)

gg_hepg2 <- gg_hepg2 +
  annotate("text", x = 0.18, y = 0.1,
           label = "Present Analysis\n(Protein-Coding Genes)",
           fontface = "bold", size = 4.5) +
  annotate("text", x = 0.82, y = 0.1,
           label = "Published\n(Ng et al., 2024)",
           fontface = "bold", size = 4.5)

# combine with patchwork
gg_islets + gg_hepg2 +
  plot_annotation(
    title = "GO Biological Process Term Concordance: Present Analysis vs. Published",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13)))