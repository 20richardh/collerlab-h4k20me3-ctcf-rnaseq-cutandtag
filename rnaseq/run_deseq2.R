# Starting point for DESeq2 differential gene expression analysis.
#
# Usage:
#   1. Edit config.R to set PROJECT_DIR and SCRIPTS_DIR.
#   2. Design DEG analysis in deg_designs/ (see DESIGN_R below).
#   3. Run.
#   4. Post-analysis: heatmap of differentially expressed genes in quiescence and 
#         chromatin remodeling GO category
#
# Outputs are written to <ANALYSIS_DIR>/<experiment_subdir>/ as defined in
# config.R.

source(file.path(dirname(sys.frame(1)$ofile), "config.R"))

library(tidyverse)

# --- Study design ------------------------------------------------------------
# Select which comparison to run by setting DESIGN_R to one of:
#   design.between_perturbations.R  (A196 vs DMSO, CTCF_OE vs RFP_OE, siCTCF vs siControl)
#   design.within_perturbation.R    (SS vs P within a single perturbation)
cutoff_name <- "sig_and_lfc"
SKIP_PCA    <- FALSE
DESIGN_R    <- "design.between_perturbations.R"

source(file.path(SCRIPTS_DIR, "deg_designs", DESIGN_R))

# --- Sample metadata ---------------------------------------------------------
DATA_COLS <- c("bioReplicate", "quiescence", "cell_type", "gene_expression", "is_OE")

select_samples <- function(df, colname, features) df[df[[colname]] %in% features, ]
exclude_samples <- function(df, colname, features) df[!(df[[colname]] %in% features), ]
select_comparison <- function(df, compare_by) {
  select_samples(df, compare_by$col, c(compare_by$ref[1], compare_by$perturb[1]))
}
index_gene_expression <- function(samples) {
  key_to_value <- GENE_EXPR_COLS
  samples %>%
    mutate(gene_expression_idx = mapply(function(x) key_to_value[[x]], gene_expression)) %>%
    mutate(is_OE = mapply(function(x) ifelse(x >= 5, "OE", "not_OE"), gene_expression_idx))
}

samples <- read.table(sample_metadata_path, header=TRUE, sep=",") %>%
  mutate(bioReplicate_id = str_replace_all(bioReplicate, "-", "_"))

if (xor(length(features_of_interest), nchar(col_of_interest)) ||
    (length(features_of_interest) > 0 &&
     !all(features_of_interest %in% samples[[col_of_interest]]))) {
  stop("DEG comparison design is improperly specified (check col_of_interest and features_of_interest)")
}

if (length(features_of_interest) > 0)
  samples <- select_samples(samples, col_of_interest, features_of_interest)
if (length(to_exclude) > 0)
  for (target in to_exclude)
    samples <- exclude_samples(samples, target$col, target$features)

samples <- samples %>%
  mutate(quiescence_idx = case_when(
    cell_type == "P"  ~ 0,
    cell_type == "CI" ~ 1,
    cell_type == "SS" ~ 2,
    TRUE ~ NA_real_
  )) %>%
  index_gene_expression()
samples$cell_type      <- factor(samples$cell_type)
samples$gene_expression <- factor(samples$gene_expression)
samples <- samples[order(samples[["quiescence_idx"]]), ]
samples <- samples[order(samples[[sort_col]]), ]
samples <- select_comparison(samples, compare_by)

# --- Build output directory --------------------------------------------------
if (length(features_of_interest) == 0) {
  subdir <- "all_experiments"
} else {
  subdir <- paste(features_of_interest, sep="_", collapse="_")
}
if (length(to_exclude) > 0) {
  excluded_str <- paste(
    sapply(to_exclude, function(t) paste(t$features, collapse="_")),
    collapse="_"
  )
  subdir <- paste(subdir, paste0("excluding_", excluded_str), sep=".")
}
if (nchar(cutoff_name)) subdir <- paste(subdir, cutoff_name, sep=".")

write_dir <- file.path(ANALYSIS_DIR, subdir)
dir.create(write_dir, showWarnings=FALSE, recursive=TRUE)
setwd(write_dir)

analysis_name <- paste(subdir, comparison_name, sep=".")

# --- Salmon quantification files --------------------------------------------
files <- samples$path
names(files) <- samples$run
stopifnot(all(file.exists(files)))

# --- Run workflow ------------------------------------------------------------
source(file.path(SCRIPTS_DIR, "deseq2_workflow.R"))
