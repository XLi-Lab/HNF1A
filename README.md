# Distinguishing Direct from Indirect HNF1A Binding

A motif-informed re-analysis of HNF1A ChIP-seq peaks and disease variants.

This repository contains the code accompanying an MSc Applied Bioinformatics
dissertation (King's College London). It re-analyses publicly available HNF1A
ChIP-seq data from Ng *et al.* (2024, *Nature Communications*; GEO **GSE206240**)
across two tissues, and asks whether the presence or absence of the HNF1A
sequence motif within a peak changes the biological interpretation drawn from it.

## Rationale

ChIP-seq identifies where a transcription factor is located, but not whether it
is bound to DNA directly through its own recognition sequence or recruited
indirectly through a protein partner. Peaks are routinely interpreted as though
all binding were direct. Here, each peak is classified by whether it contains the
HNF1A motif (JASPAR **MA0046.1**), and the two classes are carried separately
through every downstream analysis: genomic-feature enrichment, GO enrichment,
cofactor discovery, ChIP-seq signal enrichment, cross-tissue sequence comparison,
and GWAS variant prioritisation.

### Datasets

| Accession | Tissue
|---|---|---|---|---|
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

Small inputs are committed directly:

- `GSM6248576_Islets_HNF1A_ab96777_peaks.bed.gz` — islet peaks (GEO GSE206240)
- `GSM6248577_HepG2_HNF1A_ab96777_peaks.bed.gz` — HepG2 peaks (GEO GSE206240)
- `JASPAR_HNF1A_MA0046.1.meme` — the primary HNF1A matrix used throughout
- Eight partner-TF matrices used in the cofactor analysis: `HNF1B` (MA0153.2),
  `HNF4A` (MA0114.5), `HNF4G` (MA0484.3), `ONECUT1` (MA0679.3), `FOXA2`
  (MA0047.4), `FOXA3` (MA1683.2), `FOSL1` (MA0477.3), `JDP2` (MA0656.2)

**Decompress the peak files before running** — the scripts read the plain `.bed`:

Inputs *not* committed, and where to obtain them:

| Input | Source | Used by |
|---|---|---|
| Reference genome (hg38) | `BSgenome.Hsapiens.UCSC.hg38` (Bioconductor R package) | stage 03 |
| Ng *et al.* supplementary (SD2 gene lists and GO terms) | Paper supplementary information, *Nature Communications* | stage 02 |
| STRING interaction export for HNF1A | https://string-db.org (v12.0) | stage 05 |
| GWAS credible sets (Parquet) | Open Targets Platform, release 26.03 | stage 08 |
| Protein-coding gene annotation | GENCODE release 44 (`gencode.v44.annotation.gtf.gz`) | stage 08 |
| GWAS trait labels | GWAS Catalog REST API | stage 08 |

The Open Targets credible-set data is far too large to commit (several GB), so
stage 08 cannot be re-run from a clone; see the note on that stage below.

### `environment/`

- `software_versions.md` — human-readable summary of languages, packages, tools,
  and reference-data releases
- `r_sessioninformation.txt` — full `sessionInfo()` capture for the R analyses
- `conda_environment.yml` — authoritative capture of the `HNF1A_04` conda
  environment used for the GWAS step on the HPC
- `python_requirements.txt` — the GWAS notebook's direct Python dependencies

## Pipeline

The analysis alternates between R scripts run locally and manual steps performed
on the MEME Suite web server. Manual checkpoints are marked ⏸ below.

| Stage | Purpose | Key output |
|---|---|---|
| `00_data_preparation` | Builds the master one-row-per-peak table: genomic annotation (ChIPseeker), nearest gene and nearest protein-coding gene with TSS distances, genes within ±10 kb, peak score | `*_peakcentric.csv` |
| `01_genomic_feature_enrichment` | Genomic-feature distribution of peaks vs a genome background; global chi-square and per-category Fisher odds ratios | `*_Enrichment_Profile.csv` |
| `02_benchmarking` | Reproduces the published GO enrichment from Ng *et al.* as a baseline, by both a protein-coding gene list and the same gene list used in the paper | GO overlap tables and figures |
| `03_motif_scanning` | Exports peak sequences to FASTA (hg38); ⏸ **FIMO** (MA0046.1) on the web server; parses hits, collapses both-strand and overlapping hits into distinct motif sites, labels each peak motif-containing or motif-lacking | `*_FIMO_AllPeaks.csv` |
| `04_GO_by_motif` | GO BP enrichment restricted to the nearest protein-coding genes of motif-containing peaks, and its overlap with the all-peak GO | `*_ProteinCodingGO_Motif.csv` |
| `05_cofactor` | Splits peak sequences into motif-containing and motif-lacking FASTAs for ⏸ **XSTREME**; derives candidate partners from STRING; ⏸ **FIMO** for each partner TF; tests partner-motif frequency between the two peak classes (Fisher, BH-adjusted) | Cofactor enrichment table |
| `06_peak_score` | Tests ChIP-seq signal enrichment against motif presence (Mann-Whitney U) and motif count (Kruskal-Wallis, Dunn's, Spearman) | Statistics and figures |
| `07_tissue_sequence_comparison` | Compares matched motif sequences between tissues; shared vs tissue-specific target genes; Hamming distance from the consensus `GGTTAATNATTAAC` | `*_Tissue_Sequence_Comparison.csv` |
| `08_gwas` | Overlaps fine-mapped GWAS credible-set variants with peaks, assigns nearest protein-coding genes, and identifies variants falling inside an HNF1A motif | Variant tables and funnel figures |

### Script naming

- `GSM6248576_*` — pancreatic islets
- `GSM6248577_*` — HepG2
- `GSM624857[6-7]_*` — scripts operating on both tissues together (typically the
  side-by-side or top-and-bottom figure scripts)

Files ending `_figure_*` produce dissertation figures; the analysis and its
statistics live in the corresponding non-figure script.

## Running the code

### R stages (00-07)

The R scripts use paths relative to a single working directory:

```
~/Downloads/GSE206240_RAW/HNF1A/
```

with per-tissue output folders (`HNF1A_Islets_GSM6248576/GSM6248576_Outputs/`,
`HNF1A_HepG2_GSM6248577/GSM6248577_Outputs/`, and `GSM624857[6-7]_Outputs/`).
To run the pipeline elsewhere, set that working directory accordingly. Note that
`00_data_preparation` reads the raw peak BED via an absolute path, which needs
editing to match your own location.

Stages are intended to be run in numerical order within a single R session.
Most scripts read their inputs from disk and can be run independently, but two
depend on objects created earlier in the session rather than on files:

- `01_genomic_feature_enrichment` uses the peak table built by stage 00
- `03_motif_scanning/GSM624857*_FIMO.r` likewise uses the stage 00 peak table

Run stage 00 first in the same session, or read the `*_peakcentric.csv` it
writes into a variable of the same name.

### GWAS stage (08)

This stage was run on the KCL CREATE cluster under SLURM, and the paths in
`intersect.sh`, `closest.sh`, and `HNF1A_Matrix1.ipynb` are absolute
`/scratch/...` paths reflecting that environment. They are not intended to be
portable, and because the Open Targets credible-set input is several gigabytes it
is not committed here. This stage is included as a record of the analysis rather
than as a runnable pipeline.

Order within the stage:

1. `intersect.sh` — sorts the credible-set BED and intersects it with each tissue's peaks (`bedtools intersect`) → `hits_<Tissue>.bed`
2. `HNF1A_Matrix1.ipynb` — filters to variants with posterior probability ≥ 0.5 and deduplicates
3. `closest.sh` — assigns the nearest protein-coding gene from the GENCODE annotation (`bedtools closest`) → `nearest_gene_<Tissue>.bed`
4. `HNF1A_Matrix1.ipynb` — determines which variants fall inside an HNF1A motif, retrieves GWAS traits, and produces the final variant table and figures

## Data Availability

All primary data used here is publicly available. 
ChIP-seq peaks are from GEO accession [GSE206240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206240).
Motif matrices are from [JASPAR](https://jaspar.elixir.no/). 
GWAS credible sets are from the [Open Targets Platform](https://platform.opentargets.org/).

## Citation

Ng, N. H. J. *et al.* (2024) *Nature Communications*
This is the source of the re-analysed ChIP-seq data.
