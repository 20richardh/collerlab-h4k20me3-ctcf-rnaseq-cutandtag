# Reference Data

This directory contains small reference and metadata files required to run the analysis scripts.
Large files (genome FASTA, BAM files, FASTQ files, Bowtie2 index) are not included here — see the
main README for download instructions and set `PROJECT_DIR` in the relevant `config.R` / `config.sh`.

## Directory layout

```
data/
  reference/
    Q_P_deseq2_CI_SS_common_122320.csv   Consensus DEG list (P/Q quiescence comparison)
    chromatin_remodelers.tsv             Chromatin remodeling gene list
    GRCh38_unified_blacklist.bed         ENCODE GRCh38 unified blacklist regions

  rnaseq/
    samples.csv                          Sample metadata for Salmon / DESeq2
    sample_col.csv                       Column annotation table for chromatin pheatmap
    norm_log2_counts_final_noquant.csv   Normalized log2 counts (all genes, no-quant filter)
    P.CI.SS.CIR.SSR.norm_counts.csv     Triplicate normalized counts (5 conditions)
    deseq2_results/                      Pre-computed per-experiment DESeq2 result CSVs
      RFP_OE.sig_and_lfc/
        RFP_OE.sig_and_lfc.SS_vs_P.DEgenes.results.csv
      DMSO.sig_and_lfc/
        DMSO.sig_and_lfc.SS_vs_P.DEgenes.results.csv
      SiControl.sig_and_lfc/
        SiControl.sig_and_lfc.SS_vs_P.DEgenes.results.csv

  cut_and_tag/
    samplesheet.q_vs_p.csv              DiffBind samplesheet (P vs Q comparison)
```

## Sources

| File | Origin |
|------|--------|
| `Q_P_deseq2_CI_SS_common_122320.csv` | Produced by DESeq2 analysis (consensus across replicates) |
| `chromatin_remodelers.tsv` | Manually curated gene list |
| `GRCh38_unified_blacklist.bed` | ENCODE blacklist v2 (Amemiya et al. 2019) |
| `samples.csv` | Sample metadata table for this study |
| `sample_col.csv` | Column annotation for pheatmap visualization |
| `norm_log2_counts_final_noquant.csv` | Normalized counts from DESeq2 (triplicate experiment) |
| `P.CI.SS.CIR.SSR.norm_counts.csv` | Normalized counts (5-condition triplicate experiment) |
| `deseq2_results/` | DESeq2 differential expression output CSVs |
| `samplesheet.q_vs_p.csv` | DiffBind samplesheet for P vs Q CTCF/H4K20me3 comparison |
