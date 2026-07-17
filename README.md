# Distinguishing Direct from Indirect HNF1A Binding

A motif-informed re-analysis of HNF1A ChIP-seq peaks and disease variants.

MSc Applied Bioinformatics research project, King's College London, 2025-26.

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
├── data/          # committed inputs (peak BEDs, JASPAR matrices, MEME hand-offs)
├── environment/   # software and package provenance
├── scripts/       # analysis pipeline, stages 00-08, plus exploratory work
├── outputs/       # intermediates, result tables, and figures
├── reem/          # B1H cross-comparison (Reem's sequences vs. ChIP-seq peaks)
└── X/             # scratch and presentation material
```

### `data/`

Small inputs are committed directly, including the outputs of the manual MEME Suite steps.
Per-tissue inputs sit in a directory named for their accession.

```
data/
├── FIMO_Outputs/                                      # partner-TF FIMO scans (islet peaks)
│   ├── HNF1B_fimo.tsv
│   ├── HNF4A_fimo.tsv
│   ├── HNF4G_fimo.tsv
│   ├── ONECUT1_fimo.tsv
│   ├── FOXA2_fimo.tsv
│   ├── FOXA3_fimo.tsv
│   ├── FOSL1_fimo.tsv
│   └── JDP2_fimo.tsv
├── JASPAR_meme_Files/                                 # matrices uploaded to the MEME web server
│   ├── JASPAR_HNF1A_MA0046.1.meme                     # the primary matrix, used throughout
│   ├── JASPAR_HNF1B_MA0153.2.meme
│   ├── JASPAR_HNF4A_MA0114.5.meme
│   ├── JASPAR_HNF4G_MA0484.3.meme
│   ├── JASPAR_ONECUT1_MA0679.3.meme
│   ├── JASPAR_FOXA2_MA0047.4.meme
│   ├── JASPAR_FOXA3_MA1683.2.meme
│   ├── JASPAR_FOSL1_MA0477.3.meme
│   └── JASPAR_JDP2_MA0656.2.meme
├── GSM6248576/                                        # pancreatic islets
│   ├── GSM6248576_Islets_HNF1A_ab96777_peaks.bed.gz   # islet peaks (GEO GSE206240)
│   ├── GSM6248576_MA0046.1.tsv                        # islet HNF1A FIMO output
│   ├── GSM6248576_SD2_GeneList.xlsx                   # Ng et al. 2024 Supplementary Data 2
│   └── GSM6248576_SD2_GO.xlsx                         # Ng et al. 2024 published GO terms
└── GSM6248577/                                        # HepG2
    ├── GSM6248577_HepG2_HNF1A_ab96777_peaks.bed.gz
    ├── GSM6248577_MA0046.1.tsv
    ├── GSM6248577_SD2_GeneList.xlsx
    └── GSM6248577_SD2_GO.xlsx
```

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
`X_exploratory/` is not part of the run order and is described separately.

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
│   ├── GSM6248576_samegenelistGO_SD2GO_overlap.r
│   ├── GSM6248577_samegenelistGO_SD2GO_overlap.r
│   ├── GSM6248576_proteincodingGO_SD2GO_overlap.r
│   ├── GSM6248577_proteincodingGO_SD2GO_overlap.r
│   ├── GSM624857[6-7]_proteincodingGO_SD2GO_overlap_figure_sidebyside.r
│   └── GSM624857[6-7]_samegenelistGO_SD2GO_overlap_figure_sidebyside.r
├── 03_motif_scanning/
│   ├── GSM6248576_exportingpeaksequences-FASTA.r
│   ├── GSM6248577_exportingpeaksequences-FASTA.r
│   │   ⏸ FIMO (MA0046.1) Using MEME Suite Web Server
│   ├── GSM6248576_FIMO.r
│   ├── GSM6248577_FIMO.r
│   └── GSM624857[6-7]_motifpresence_figure_sidebyside.r
├── 04_GO_by_motif/
│   ├── GSM6248576_motifcontainingGO.r
│   ├── GSM6248576_proteincodingGO_proteincodingmotifGO_overlap.r
│   └── GSM6248576_proteincodingGO_proteincodingmotifGO_overlap_figure_sidebyside.r
├── 05_cofactor/
│   ├── GSM6248576_HNF1A_with_withoutMOTIF.r
│   │   ⏸ XSTREME on Motif-Lacking FASTA Using MEME Suite Web Server
│   ├── HNF1A_STRING_interactors.r
│   │   ⏸ FIMO For Each Partner-TF Matrix Using MEME Suite Web Server
│   └── GSM6248576_interactionpartners_enrichment.r
├── 06_peak_score/
│   ├── GSM6248576_peakscore.r
│   ├── GSM6248577_peakscore.r
│   ├── GSM624857[6-7]_peakscore_binary_figure_sidebyside_significance.r
│   └── GSM624857[6-7]_peakscore_counts_figure_sidebyside_significance.r
├── 07_tissue_sequence_comparison/
│   ├── GSM624857[6-7]_FIMO_comparison.r
│   └── GSM624857[6-7]_PPM_comparison_figure_topandbottom.r
├── 08_gwas/
│   ├── HNF1A_Matrix1.ipynb
│   ├── intersect.sh
│   └── closest.sh
└── X_exploratory/
    ├── GSM6248576/
    │   ├── GSM6248576_FIMO_motifpeakcounts.r
    │   ├── GSM6248576_GOanalysis.r
    │   └── GSM6248576_motifGO_SD2GO_overlap.r
    ├── GSM6248577/
    │   ├── GSM6248577_GOanalysis.r
    │   ├── GSM6248577_HNF1A_with_withoutMOTIF.r
    │   └── GSM6248577_proteincodingGO_proteincodingmotifGO_overlap.r
    └── GSM624857[6-7]/
        └── GSM624857[6-7]_coordvsgene_specificity.r
```

Stages `04` and `05` are islet-only in the run order. 
The HepG2 counterparts were run but are held in `X_exploratory/` rather than the pipeline.

### `outputs/`

Intermediates and results, in directories named for the accession they derive from. 

```
outputs/
├── GSM6248576_Outputs/                          # pancreatic islets
│   ├── GSM6248576_PeakCentric.csv               # stage 00
│   ├── GSM6248576_Enrichment_Profile.csv        # stage 01
│   ├── GSM6248576_SameGeneListGO.csv            # stage 02
│   ├── GSM6248576_ProteinCodingGO_AllPeaks.csv  # stage 02
│   ├── GSM6248576_peaks-FASTAs.fa               # stage 03, input to FIMO
│   ├── GSM6248576_FIMO_AllPeaks.csv             # stage 03
│   ├── GSM6248576_ProteinCodingGO_Motif.csv     # stage 04
│   ├── GSM6248576_HNF1A_withMotif.fa            # stage 05
│   ├── GSM6248576_HNF1A_withoutMotif.fa         # stage 05, input to XSTREME
│   ├── GSM6248576_SD2_GeneList.xlsx             # copy of the data/ input
│   └── GSM6248576_SD2_GO.xlsx                   # copy of the data/ input
├── GSM6248577_Outputs/                          # HepG2
│   ├── GSM6248577_PeakCentric.csv
│   ├── GSM6248577_Enrichment_Profile.csv
│   ├── GSM6248577_SameGeneListGO.csv
│   ├── GSM6248577_ProteinCodingGO_AllPeaks.csv
│   ├── GSM6248577_peaks-FASTAs.fa
│   └── GSM6248577_FIMO_AllPeaks.csv
├── GSM624857[6-7]_Outputs/
│   └── GSM624857[6-7]_Tissue_Sequence_Comparison.csv   # stage 07
├── Figures/
│   ├── JASPAR_MATRIX.png
│   ├── genomicenrichment.png
│   ├── allvsmotifcomparison_0.02cutoff.png
│   ├── motifpresence.png
│   ├── motifcount.png
│   ├── isletsmotif1.png
│   ├── hepg2motif2.png
│   ├── motifpresenceisletsbinary.png
│   ├── motifpresenceisletssequential.png
│   ├── motifpresencehepbinary.png
│   ├── motifpresencehepsequential.png
│   ├── motifpresenceisletshep.png
│   └── logohepislets1.png
└── X/
    ├── HNF1A_STRING_interactors.tsv             # stage 05
    ├── HNF1A_Islets_common.xlsx
    └── HNF1A_HepG2_common.xlsx
```

The HepG2 output set is smaller than the islet set because stages `04` and `05` are islet-only.

### `reem/`

Cross-comparison between Reem's bacterial one-hybrid (B1H) sequences and the ChIP-seq peaks. 
Each B1H sequence is searched against the peak FASTAs from stage `03` on both strands, and matched peaks are joined to their ChIP-seq peak scores. 
Two match stringencies are used throughout: the exact 20 bp sequence, and the 14 bp core with the ±3 bp flanks removed.

```
reem/
├── reem_enrichmentscore_peakscore_1.r   # ΔDBD − Lib vs. peak score
├── reem_enrichmentscore_peakscore_2.r
├── reem_DBD_log2FC_mean.r               # DBD_log2FC_mean vs. peak score
├── reem_Lib_log2FC_mean.r               # Lib_log2FC_mean vs. peak score
├── reem_delta_DBD_minus_Lib.r
├── reem1.R                              # scratch
├── CSVs/
│   ├── HNF1A_B1H_all_motifs_corrected_2.csv   # Reem's B1H table (source input)
│   ├── reem_peak_matches.csv
│   ├── reem_exact20bp_matches.csv
│   ├── reem_core14bp_matches.csv
│   ├── reem_DBD_exact20bp_matches.csv
│   ├── reem_DBD_core14bp_matches.csv
│   ├── reem_Lib_exact20bp_matches.csv
│   └── reem_Lib_core14bp_matches.csv
└── Plots/
    ├── reem_peak_scatter.pdf
    ├── reem_scatter_exact20bp.pdf
    ├── reem_scatter_core14bp.pdf
    ├── reem_DBD_scatter_exact20bp.pdf
    ├── reem_DBD_scatter_core14bp.pdf
    ├── reem_Lib_scatter_exact20bp.pdf
    └── reem_Lib_scatter_core14bp.pdf
```

### `X/`

Scratch and presentation material. Not part of the analysis.


## Pipeline

The analysis alternates between R scripts run locally and manual steps performed on the MEME Suite web server.

| Stage |
|---|
| `00_data_preparation` |
| `01_genomic_feature_enrichment` |
| `02_benchmarking` |
| `03_motif_scanning` |
| `04_GO_by_motif` |
| `05_cofactor` |
| `06_peak_score` |
| `07_tissue_sequence_comparison` |
| `08_gwas` |

## Citation

The re-analysed ChIP-seq data are from:

Ng, N. H. J., Ghosh, S., Bok, C. M., Ching, C., Low, B. S. J., Chen, J. T., Lim, E.,
Miserendino, M. C., Tan, Y. S., Hoon, S., & Teo, A. K. K. (2024). HNF4A and HNF1A
exhibit tissue specific target gene regulation in pancreatic beta cells and
hepatocytes. *Nature Communications*, 15(1), 4288.
https://doi.org/10.1038/s41467-024-48647-w
