# =============================================================================
# USER CONFIGURATION — edit this file before running any RNA-seq scripts.
#
# Directory layout assumed under PROJECT_DIR:
#   genome/                  genome annotation files (Gencode v29, hg38)
#   salmon_quant/            Salmon output folders (one per sample)
#   salmon_quant/samples.csv sample metadata table
#   deseq/                   reference tables (chromatin remodeler list, etc.)
#   results/DESeq2/          analysis outputs will be written here
#
# Reference files produced by a prior analysis (triplicate experiment) that
# are read by heatmaps.R and chromatin_pheatmap.R must also live under
# PROJECT_DIR (see paths below).
# =============================================================================

PROJECT_DIR <- "/path/to/your/project"         # root data directory
SCRIPTS_DIR <- "/path/to/this/repository/rnaseq"  # path to THIS scripts folder

# --- Derived paths (adjust only if your layout differs) ----------------------
GENOME_DIR    <- file.path(PROJECT_DIR, "genome")
SALMON_DIR    <- file.path(PROJECT_DIR, "salmon_quant")
ANALYSIS_DIR  <- file.path(PROJECT_DIR, "results", "DESeq2")
DESEQ_DIR     <- file.path(PROJECT_DIR, "deseq")

# Genome annotation (Gencode v29 / hg38)
genome_path      <- file.path(GENOME_DIR, "gencode.v29.annotation.gff3")
genome_csv_path  <- file.path(GENOME_DIR, "gencode.v29.annotation.csv")

# Sample metadata CSV (columns: run, path, cell_type, gene_expression, bioReplicate, ...)
sample_metadata_path <- file.path(SALMON_DIR, "samples.csv")

# Reference data files
CONSENSUS_DEG_FILE        <- file.path(PROJECT_DIR, "Q_P_deseq2_CI_SS_common_122320.csv")
CHROMATIN_REMODELERS_FILE <- file.path(DESEQ_DIR,   "chromatin_remodelers.tsv")
NORM_COUNTS_FILE          <- file.path(DESEQ_DIR,   "norm_log2_counts_final_noquant.csv")
TRIPLICATE_NORM_COUNTS_FILE <- file.path(DESEQ_DIR, "P.CI.SS.CIR.SSR.norm_counts.csv")

# Per-experiment DESeq2 result CSVs used by heatmaps.R
DESEQ_RESULTS_DIR <- file.path(PROJECT_DIR, "results", "DESeq2")

# --- Analysis parameters -----------------------------------------------------
THRESHOLD_P   <- 0.05
THRESHOLD_LFC <- 1

# Column ordering for gene expression conditions (used for PCA coloring)
GENE_EXPR_COLS <- list(DMSO=1, A196=2, SiControl=3, SiCTCF=4, RFP_OE=5, CTCF_OE=6)
