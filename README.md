# Distinguishing Direct from Indirect HNF1A Binding

A motif-informed re-analysis of HNF1A ChIP-seq peaks and disease variants.

## Data Availability

All primary data used here is publicly available. 
ChIP-seq peaks are from GEO accession [GSE206240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206240).
Motif matrices are from [JASPAR](https://jaspar.elixir.no/). 
GWAS credible sets are from the [Open Targets Platform](https://platform.opentargets.org/).

### Datasets

| Accession | Tissue
|---|---|---|
| GSM6248576 | Primary human pancreatic islets
| GSM6248577 | HepG2 hepatocellular carcinoma

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
- Eight partner-TF matrices used in the cofactor analysis: `HNF1B` (MA0153.2),
  `HNF4A` (MA0114.5), `HNF4G` (MA0484.3), `ONECUT1` (MA0679.3), `FOXA2`
  (MA0047.4), `FOXA3` (MA1683.2), `FOSL1` (MA0477.3), `JDP2` (MA0656.2)

Inputs *not* committed due to size, and where to obtain them.

| Input | Source | Used By |
|---|---|---|
| GWAS Credible Sets (Parquet) | Open Targets Platform, release 26.03 | stage 08 |
| GWAS Trait Labels | GWAS Catalog REST API | stage 08 |

### `environment/`

- `software_versions.md`
- `r_sessioninformation.txt`
- `conda_environment.yml`
- `python_requirements.txt`

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

### Script Naming

- `GSM6248576_*` — pancreatic islets
- `GSM6248577_*` — HepG2
- `GSM624857[6-7]_*` — scripts operating on both tissues together

Files ending `_figure_*` produce dissertation figures.

## Citation

Ng, N. H. J. *et al.* (2024) *Nature Communications*
This is the source of the re-analysed ChIP-seq data.
