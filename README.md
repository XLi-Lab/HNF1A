# Distinguishing Direct from Indirect HNF1A Binding

A motif-informed re-analysis of HNF1A ChIP-seq peaks and disease variants.

## Data Availability

All primary data used here is publicly available. 
ChIP-seq peaks are from GEO accession [GSE206240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206240).
Motif matrices are from [JASPAR](https://jaspar.elixir.no/). 
GWAS credible sets are from the [Open Targets Platform](https://platform.opentargets.org/).

### Datasets

| Accession | Tissue
|---|---|
| GSM6248576 | Primary Human Pancreatic Islets
| GSM6248577 | HepG2 Hepatocellular Carcinoma

## Repository Structure

```
.
├── data/          # committed inputs (peak BEDs, JASPAR motif matrices)
├── environment/   # software and package provenance
└── scripts/       # analysis pipeline, stages 00-08
```

### `data/`

Small inputs are committed directly.

- `GSM6248576_Islets_HNF1A_ab96777_peaks.bed.gz` — islet peaks (GEO GSE206240)
- `GSM6248577_HepG2_HNF1A_ab96777_peaks.bed.gz` — HepG2 peaks (GEO GSE206240)
- `JASPAR_HNF1A_MA0046.1.meme` — the primary HNF1A matrix used throughout
- Eight partner-TF matrices used in the cofactor analysis — 
`HNF1B` (MA0153.2),
`HNF4A` (MA0114.5),
`HNF4G` (MA0484.3),
`ONECUT1` (MA0679.3),
`FOXA2` (MA0047.4),
`FOXA3` (MA1683.2),
`FOSL1` (MA0477.3),
`JDP2` (MA0656.2)

Inputs *not* committed due to size, and where to obtain them.

| Input | Source | Used By |
|---|---|---|
| GWAS Credible Sets (Parquet) | Open Targets Platform, Release 26.03 | Stage 08 |
| GWAS Trait Labels | GWAS Catalog REST API | Stage 08 |

### `environment/`

- `software_versions.md`
- `r_sessioninformation.txt`
- `conda_environment.yml`
- `python_requirements.txt`

### `scripts/`

Script names carry the GEO accession of the tissue they operate on.
- `GSM6248576_*` — pancreatic islets
- `GSM6248577_*` — HepG2
- `GSM624857[6-7]_*` — both tissues together
Files ending `_figure_*` produce dissertation figures.
Scripts are listed below in run order.

```
scripts/
├── 00_data_preparation/
│   ├── GSM6248576_peakcentric_updated.r
│   └── GSM6248577_peakcentric_updated.r
├── 01_genomic_feature_enrichment/
│   ├── GSM6248576_enrichmentprofile.r
│   ├── GSM6248577_enrichmentprofile.r
│   └── GSM624857[6-7]_enrichmentprofile_figure_topandbottom.r
├── 02_benchmarking/
│   ├── GSM6248576_proteincodingGO_SD2GO_overlap.r
│   ├── GSM6248577_proteincodingGO_SD2GO_overlap.r
│   ├── GSM6248576_samegenelistGO_SD2GO_overlap.r
│   ├── GSM6248577_samegenelistGO_SD2GO_overlap.r
│   ├── GSM624857[6-7]_proteincodingGO_SD2GO_overlap_figure_sidebyside.r
│   └── GSM624857[6-7]_samegenelistGO_SD2GO_overlap_figure_sidebyside.r
├── 03_motif_scanning/
│   ├── GSM6248576_exportingpeaksequences-FASTA.r
│   ├── GSM6248577_exportingpeaksequences-FASTA.r
│   │   ⏸ FIMO (MA0046.1) on the MEME Suite web server
│   ├── GSM6248576_FIMO.r
│   ├── GSM6248577_FIMO.r
│   └── GSM624857[6-7]_motifpresence_figure_sidebyside.r
├── 04_GO_by_motif/
│   ├── GSM6248576_motifcontainingGO.r
│   ├── GSM6248576_proteincodingGO_proteincodingmotifGO_overlap.r
│   └── GSM6248576_proteincodingGO_proteincodingmotifGO_overlap_figure_sidebyside.r
├── 05_cofactor/
│   ├── GSM6248576_HNF1A_with_withoutMOTIF.r
│   │   ⏸ XSTREME on the motif-lacking FASTA (web server)
│   ├── HNF1A_STRING_interactors.r
│   │   ⏸ FIMO for each partner-TF matrix (web server)
│   └── GSM6248576_interactionpartners_enrichment.r
├── 06_peak_score/
│   ├── GSM6248576_peakscore.r
│   ├── GSM6248577_peakscore.r
│   ├── GSM624857[6-7]_peakscore_binary_figure_sidebyside_significance.r
│   └── GSM624857[6-7]_peakscore_counts_figure_sidebyside_significance.r
├── 07_tissue_sequence_comparison/
│   ├── GSM624857[6-7]_FIMO_comparison.r
│   └── GSM624857[6-7]_PPM_comparison_figure_topandbottom.r
└── 08_gwas/
    ├── intersect.sh
    ├── HNF1A_Matrix1.ipynb
    └── closest.sh
```

## Pipeline

The analysis alternates between R scripts run locally and manual steps performed on the MEME Suite web server.

| Stage | Key Output |
|---|---|
| `00_data_preparation` | `*_peakcentric.csv` |
| `01_genomic_feature_enrichment` | `*_Enrichment_Profile.csv` |
| `02_benchmarking` | GO Overlap Tables & Figures |
| `03_motif_scanning` | `*_FIMO_AllPeaks.csv` |
| `04_GO_by_motif` | `*_ProteinCodingGO_Motif.csv` |
| `05_cofactor` | Cofactor Enrichment Table |
| `06_peak_score` | Statistics & Figures |
| `07_tissue_sequence_comparison` | `*_Tissue_Sequence_Comparison.csv` |
| `08_gwas` | Variant Tables & Figures |

## Citation

Ng, N. H. J. *et al.* (2024) *Nature Communications*
This is the source of the re-analysed ChIP-seq data.
