# =============================================================================
# USER CONFIGURATION — edit this file before running diffbind_workflow.R.
#
# Set PROJECT_DIR to your large data directory:
#   peakCalling/SEACR/   SEACR peak files (*_seacr_NoBg_T10.peaks.stringent.no_blacklist.23chr_only.bed)
#   bam/CTCF/            BAM files for CTCF
#   bam/H4K20me3/        BAM files for H4K20me3
#   genome/              gencode.v29.annotation.csv (Gencode v29)
#
# The DiffBind samplesheet and genome blacklist are bundled under data/ in
# this repository and are auto-detected — no path changes needed for those.
# =============================================================================

PROJECT_DIR <- "/path/to/your/project/data"

# --- Auto-detected repo paths (no change needed) -----------------------------
DIFFBIND_DIR <- dirname(sys.frame(1)$ofile)               # cut_and_tag/diffbind/
REPO_DIR     <- dirname(dirname(DIFFBIND_DIR))             # repository root
DATA_DIR     <- file.path(REPO_DIR, "data")

# --- Large external data (derived from PROJECT_DIR) --------------------------
GENOME_DIR   <- file.path(PROJECT_DIR, "genome")
GENOME_CSV   <- file.path(GENOME_DIR, "gencode.v29.annotation.csv")
BAM_DIR      <- file.path(PROJECT_DIR, "bam")
BIGWIG_DIR   <- BAM_DIR

PEAK_CALLING_DIR <- file.path(PROJECT_DIR, "peakCalling", "SEACR")
HIC_DIR          <- file.path(PROJECT_DIR, "hic")

# --- Bundled reference files (data/ in this repo) ----------------------------
CHIPSEQ_BLACKLIST <- file.path(DATA_DIR, "reference", "GRCh38_unified_blacklist.bed")
SAMPLESHEET_QP    <- file.path(DATA_DIR, "cut_and_tag", "samplesheet.q_vs_p.csv")

# --- Analysis parameters -----------------------------------------------------
BIN_SIZE   <- 200
PEAK_TYPE  <- "T10"
FDR_THRESHOLD <- 0.05

PEAK_FILE_SUFFIX <- paste0("_seacr_NoBg_", PEAK_TYPE,
                            ".peaks.stringent.no_blacklist.23chr_only.bed")
BIGWIG_SUFFIX    <- ".rmDup.mapped_only.no_blacklist.bw"

# --- Experiment metadata -----------------------------------------------------
FEATURES     <- c("CTCF", "H4K20me3")
BIOREPLICATES <- c("12-3", "12-4")
CELL_TYPES   <- c("P", "CI", "SS")
SEP2         <- "-"

CHROMOSOMES <- c(paste0("chr", 1:22), "chrX", "chrY")
