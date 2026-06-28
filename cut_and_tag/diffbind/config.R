# =============================================================================
# USER CONFIGURATION — edit this file before running diffbind_workflow.R.
#
# Expected directory layout under PROJECT_DIR:
#   peakCalling/SEACR/   SEACR peak files (*_seacr_NoBg_T10.peaks.stringent.no_blacklist.23chr_only.bed)
#   bam/CTCF/            BAM files for CTCF
#   bam/H4K20me3/        BAM files for H4K20me3
#   diffbind/metadata/   samplesheet CSVs for DiffBind
#   genome/              gencode.v29.annotation.csv
#   bed/reference/blacklist/  GRCh38_unified_blacklist.bed
# =============================================================================

PROJECT_DIR <- "/path/to/your/project/data"

# --- Derived paths -----------------------------------------------------------
GENOME_DIR    <- file.path(PROJECT_DIR, "genome")
GENOME_CSV    <- file.path(GENOME_DIR, "gencode.v29.annotation.csv")
BAM_DIR       <- file.path(PROJECT_DIR, "bam")
BIGWIG_DIR    <- BAM_DIR

CHIPSEQ_BLACKLIST <- file.path(PROJECT_DIR, "bed", "reference", "blacklist",
                                "GRCh38_unified_blacklist.bed")
PEAK_CALLING_DIR  <- file.path(PROJECT_DIR, "peakCalling", "SEACR")
HIC_DIR           <- file.path(PROJECT_DIR, "hic")

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
