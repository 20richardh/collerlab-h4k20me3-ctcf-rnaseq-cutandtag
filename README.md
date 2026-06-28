# Analysis Scripts — H4K20me3 and CTCF in Quiescence

Analysis accompanying:

> **H4K20me3 and CTCF act reciprocally at TAD boundaries to regulate cell state transitions**  
> Hilary Coller lab, UCLA  
> https://doi.org/10.64898/2026.01.29.702647

This repository contains scripts for:
1. **RNA-seq differential gene expression** (DESeq2)
2. **CUT&Tag peak calling** (SEACR pipeline) — upstream of all CUT&Tag analyses
3. **Differential CUT&Tag binding analysis** (DiffBind)


---

## Repository structure

```
data/                         # bundled reference and metadata files
  reference/
    Q_P_deseq2_CI_SS_common_122320.csv   consensus DEG list
    chromatin_remodelers.tsv             chromatin remodeling gene list
    GRCh38_unified_blacklist.bed         ENCODE GRCh38 unified blacklist
  rnaseq/
    samples.csv                          RNA-seq sample metadata
    sample_col.csv                       column annotation for pheatmap
    norm_log2_counts_final_noquant.csv   normalized log2 counts
    P.CI.SS.CIR.SSR.norm_counts.csv     triplicate normalized counts
    deseq2_results/                      pre-computed DESeq2 result CSVs
  cut_and_tag/
    samplesheet.q_vs_p.csv              DiffBind samplesheet (P vs Q)

rnaseq/
  config.R                  # set PROJECT_DIR here (genome + Salmon paths)
  run_deseq2.R              # entry point: sources config.R + a design file, runs deseq2_workflow.R
  deseq2_workflow.R         # tximport → DESeq2 → PCA → heatmap → GO enrichment
  heatmaps.R                # comparison heatmaps (Figs 1A, 3E); run in same session as run_deseq2.R
  chromatin_pheatmap.R      # chromatin remodeling gene heatmap (Fig 1B); largely standalone
  deg_designs/              # experimental design configs (one per comparison)
    design.between_perturbations.R
    design.within_perturbation.R

cut_and_tag/
  peak_calling/             # SEACR CUT&Tag peak-calling pipeline (steps 1–4, SGE cluster)
    config.sh               # set PROJ_PATH here
    step1_align.sh          # Bowtie2 alignment
    step2_dedup.sh          # Picard duplicate removal
    step3_filter.sh         # BAM filtering, BED conversion, fragment binning
    step4_peak_call.sh      # Bedgraph generation + SEACR peak calling
    SEACR_1.3.sh            # SEACR tool (Meers et al. 2019)
  diffbind/                 # Differential binding analysis (Figs 2C, 2D)
    config.R                # set PROJECT_DIR here (BAM + peak file paths)
    diffbind_workflow.R     # DiffBind + DESeq2 differential peak analysis
```

---

## Requirements

### R packages

Install from CRAN:
```r
install.packages(c("tidyverse", "ggplot2", "reshape2", "pheatmap",
                   "RColorBrewer", "ggfortify", "grid", "gridExtra"))
```

Install from Bioconductor:
```r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2", "tximport", "AnnotationDbi", "GenomicFeatures", "txdbmaker",
  "apeglm", "org.Hs.eg.db", "clusterProfiler", "PCAtools",
  "ComplexHeatmap", "InteractiveComplexHeatmap", "circlize",
  "PoiClaClu", "ChIPpeakAnno", "DiffBind", "GenomicRanges",
  "biomaRt", "fgsea"
))
install.packages("viridis")  # also on CRAN
```

### Command-line tools (peak calling pipeline)

- **Bowtie2** ≥ 2.4.2
- **SAMtools** ≥ 1.15
- **Picard Tools**
- **BEDTools** ≥ 2.30.0
- **SEACR** 1.3 (included as `SEACR_1.3.sh`; requires R)
- **deeptools** (via conda/pip, used for bigwig generation)

Genome: **hg38 / GRCh38**, annotation **Gencode v29**.

---

## Data

Raw sequencing data are deposited at GEO: **[accession number]**.

### Bundled in this repository (`data/`)

These files are included and require no additional setup:

| File | Used by |
|------|---------|
| `data/reference/Q_P_deseq2_CI_SS_common_122320.csv` | `deseq2_workflow.R`, `heatmaps.R`, `chromatin_pheatmap.R` |
| `data/reference/chromatin_remodelers.tsv` | `chromatin_pheatmap.R`, `deseq2_workflow.R` |
| `data/reference/GRCh38_unified_blacklist.bed` | `step3_filter.sh`, `step4_peak_call.sh`, `diffbind_workflow.R` |
| `data/rnaseq/samples.csv` | `run_deseq2.R` |
| `data/rnaseq/sample_col.csv` | `chromatin_pheatmap.R` |
| `data/rnaseq/norm_log2_counts_final_noquant.csv` | `chromatin_pheatmap.R` |
| `data/rnaseq/P.CI.SS.CIR.SSR.norm_counts.csv` | `heatmaps.R` |
| `data/rnaseq/deseq2_results/` (3 CSV files) | `heatmaps.R` |
| `data/cut_and_tag/samplesheet.q_vs_p.csv` | `diffbind_workflow.R` |

### Large files — set paths in config

These must be obtained externally and configured via `PROJECT_DIR` in `config.R` / `PROJ_PATH` in `config.sh`:

| File | Source | Used by |
|------|--------|---------|
| Salmon `quant.sf` files (one per sample) | GEO (this study) | `run_deseq2.R` |
| `gencode.v29.annotation.gff3` | [Gencode release 29](https://www.gencodegenes.org/human/release_29.html) | `run_deseq2.R` |
| `gencode.v29.annotation.csv` (derived) | Generated from GFF3 by workflow | `chromatin_pheatmap.R` |
| `GRCh38.primary_assembly.genome.fa` | [Gencode release 29](https://www.gencodegenes.org/human/release_29.html) | CUT&Tag pipeline |
| hg38 Bowtie2 index | Built from genome FASTA | `step1_align.sh` |
| SEACR peak BED files | Produced by `step4_peak_call.sh` | `diffbind_workflow.R` |
| BAM files (CTCF, H4K20me3) | Produced by `step2_dedup.sh` | `diffbind_workflow.R` |
| Raw FASTQ files | GEO (this study) | `step1_align.sh` |

---

## Usage

### RNA-seq (DESeq2)

1. Edit `rnaseq/config.R` — set `PROJECT_DIR` (genome + Salmon quant paths). `SCRIPTS_DIR` and reference data paths are auto-detected.
2. Choose a design in `rnaseq/deg_designs/` and set `DESIGN_R` in `run_deseq2.R`.
3. Source `run_deseq2.R` in an R session.
4. Optionally, in the **same R session** (objects still loaded), source `heatmaps.R` (Figs 1A, 3E).
5. `chromatin_pheatmap.R` (Fig 1B) can be run independently after outputs are on disk.

The design used for the primary comparisons in the paper is `design.between_perturbations.R`.

### CUT&Tag peak calling (cluster)

Run steps in order. Each is an independent SGE job:

```bash
qsub step1_align.sh
# After step 1 completes:
qsub step2_dedup.sh
qsub step3_filter.sh
qsub step4_peak_call.sh
```

Edit `config.sh` first — set `PROJ_PATH` to your project root and `EMAIL` for job notifications.
Adjust `module load` lines for your cluster's module system.

### DiffBind (Figs 2C, 2D)

1. Edit `cut_and_tag/diffbind/config.R` — set `PROJECT_DIR`.
2. Source `diffbind_workflow.R` in R.

---


## Key dependencies

- **DESeq2:** Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550. https://doi.org/10.1186/s13059-014-0550-8

- **tximport:** Soneson C, Love MI, Robinson MD (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research* 4:1521. https://doi.org/10.12688/f1000research.7563.2

- **SEACR:** Meers MP, Tenenbaum D, Henikoff S (2019). Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling. *Epigenetics & Chromatin* 12:42. https://doi.org/10.1186/s13072-019-0287-4

- **DiffBind:** Ross-Innes CS, Stark R, Teschendorff AE, Holmes KA, Ali HR, Dunning MJ, Brown GD, Gojis O, Ellis IO, Green AR, Ali S, Chin S-F, Palmieri C, Caldas C, Carroll JS (2012). Differential oestrogen receptor binding is associated with clinical outcome in breast cancer. *Nature* 481:389–393. https://doi.org/10.1038/nature10730

---

## License

MIT License. See `LICENSE`.
