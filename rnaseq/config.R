# =============================================================================
# USER CONFIGURATION — edit this file before running any RNA-seq scripts.
#
# Two roots to set:
#   PROJECT_DIR  — your large data directory (genome FASTA, Salmon quant output)
#   REPO_DIR     — auto-detected from this file's location (no change needed)
#
# Small reference files (consensus DEG list, normalized counts, samplesheets)
# are bundled under data/ in this repository and require no path changes.
#
# Large files that must be supplied externally (set PROJECT_DIR):
#   genome/gencode.v29.annotation.gff3    Gencode v29 GFF3
#   genome/gencode.v29.annotation.csv     Gencode v29 gene table (derived)
#   salmon_quant/<sample>/quant.sf        Salmon quantification output
#   results/DESeq2/                       analysis outputs are written here
# =============================================================================

PROJECT_DIR <- "/path/to/your/project"   # root for large external data

# --- Auto-detected repo paths (no change needed) -----------------------------
SCRIPTS_DIR <- dirname(sys.frame(1)$ofile)          # rnaseq/ directory
REPO_DIR    <- dirname(SCRIPTS_DIR)                  # repository root
DATA_DIR    <- file.path(REPO_DIR, "data")           # bundled reference data

# --- Large external data (derived from PROJECT_DIR) --------------------------
GENOME_DIR   <- file.path(PROJECT_DIR, "genome")
SALMON_DIR   <- file.path(PROJECT_DIR, "salmon_quant")
ANALYSIS_DIR <- file.path(PROJECT_DIR, "results", "DESeq2")  # outputs written here

# Genome annotation (Gencode v29 / hg38) — download from Gencode release 29
genome_path     <- file.path(GENOME_DIR, "gencode.v29.annotation.gff3")
genome_csv_path <- file.path(GENOME_DIR, "gencode.v29.annotation.csv")

# --- Bundled reference files (data/ in this repo) ----------------------------
# Sample metadata CSV (columns: run, path, cell_type, gene_expression, bioReplicate, ...)
sample_metadata_path <- file.path(DATA_DIR, "rnaseq", "samples.csv")

CONSENSUS_DEG_FILE          <- file.path(DATA_DIR, "reference", "Q_P_deseq2_CI_SS_common_122320.csv")
CHROMATIN_REMODELERS_FILE   <- file.path(DATA_DIR, "reference", "chromatin_remodelers.tsv")
NORM_COUNTS_FILE            <- file.path(DATA_DIR, "rnaseq",    "norm_log2_counts_final_noquant.csv")
TRIPLICATE_NORM_COUNTS_FILE <- file.path(DATA_DIR, "rnaseq",    "P.CI.SS.CIR.SSR.norm_counts.csv")

# Pre-computed per-experiment DESeq2 result CSVs read by heatmaps.R
# (also re-generated under ANALYSIS_DIR when run_deseq2.R is run from scratch)
DESEQ_RESULTS_DIR <- file.path(DATA_DIR, "rnaseq", "deseq2_results")

# sample_col.csv annotation table used by chromatin_pheatmap.R via setwd(DESEQ_DIR)
DESEQ_DIR <- file.path(DATA_DIR, "rnaseq")

# --- Analysis parameters -----------------------------------------------------
THRESHOLD_P   <- 0.05
THRESHOLD_LFC <- 1

# Column ordering for gene expression conditions (used for PCA coloring)
GENE_EXPR_COLS <- list(DMSO=1, A196=2, SiControl=3, SiCTCF=4, RFP_OE=5, CTCF_OE=6)
