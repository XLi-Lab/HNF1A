# software environment

## languages
- **R** v4.5.1 — run locally on macOS Sonoma 14.6 (peak annotation, benchmarking,
  motif-based classification, GO enrichment, peak-score, sequence analyses).
- **python** v3.10.20 — run on the King's College London CREATE HPC cluster
  (GWAS credible-set overlap), in the conda-managed environment `HNF1A_04`;
  jobs submitted via SLURM.

## key R packages
- ChIPseeker 1.44.0
- clusterProfiler 4.16.0
- GenomicRanges, IRanges, Biostrings
- BSgenome.Hsapiens.UCSC.hg38 1.4.5
- TxDb.Hsapiens.UCSC.hg38.knownGene 3.21.0
- EnsDb.Hsapiens.v86
- org.Hs.eg.db
- readxl 1.4.5 (reads the Ng et al. 2024 supplementary tables in benchmarking)
- ggplot2, patchwork, ggseqlogo 0.2.2, gridExtra, VennDiagram
- FSA (Dunn's test)
- readr, dplyr, tidyr, purrr

The full, exact package set is captured in `r_sessioninformation.txt`.

## python packages (HPC)
run in the `HNF1A_04` conda environment on CREATE (Python 3.10.20):
- pandas 2.3.3
- numpy 2.2.5
- pyarrow 23.0.1
- matplotlib 3.10.9
- requests 2.33.1

`python_requirements.txt` lists the direct dependencies of the GWAS notebook.
`conda_environment.yml` is the full, authoritative capture of the `HNF1A_04`
environment (conda + pip), exported with `conda env export --no-builds`.

## command-line/system tools
- bedtools 2.31.0 (HPC)

## MEME Suite (web server, v5.5.9)
Run interactively on the MEME Suite web server.
- **FIMO** — motif scanning of peak sequences.
- **XSTREME** — de novo motif discovery on motif-lacking islet peaks (cofactor arm).

