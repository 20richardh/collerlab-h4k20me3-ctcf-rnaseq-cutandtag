# Comparison heatmaps for Figures 1A and 3E.
# Run this script in the same R session after run_deseq2.R has completed,
# OR load a saved workspace: load("<analysis_dir>/workspace.RData")
#
# Requires workspace objects: dds, geneID_to_geneName, sampleTable,
#   sigGenes.consensus, upGenes.consensus, downGenes.consensus,
#   counts.lognorm, compare_by, THRESHOLD_P, THRESHOLD_LFC
# Also requires config.R variables for file paths.
#-------------------------------------------------------------------------------
ANALYSIS <- "A196" # OE # A196 # siCTCF
GENE_UNIVERSE <- "deg_reference" # current: taking intersect of deg_reference x triplicate genes # "deg_reference" # "triplicate"
GENE_SUBSET <- "up_in_q"
#-------------------------------------------------------------------------------

if (ANALYSIS == "OE") {
  # OE, P/SS
  COLS_MASK <- c("12-3_CI_RFP_OE_S33", "12-4_CI_RFP_OE_S34",
                 "12-3_CI_CTCF_OE_S35", "12-4_CI_CTCF_OE_S36")
  BASE_COLS <- c("12-3_P_RFP_OE_S25", "12-4_P_RFP_OE_S26")
  DIFF_COLS <- c("12-3_SS_RFP_OE_S29", "12-4_SS_RFP_OE_S30")
  BASE2 <- c("12-3_P_CTCF_OE_S35", "12-4_P_CTCF_OE_S36")
  DIFF2 <- c("12-3_SS_CTCF_OE_S31", "12-4_SS_CTCF_OE_S32")
} else if (ANALYSIS == "A196") {
  # A196, SS
  COLS_MASK <- c("12-3_CI_DMSO_S9", "12-4_CI_DMSO_S10",
                 "12-3_CI_A196_S11", "12-4_CI_A196_S12")#,
                 #"12-3_P_DMSO_S1", "12-4_P_DMSO_S2",
                 #"12-3_P_A196_S3", "12-4_P_A196_S4")
  BASE_COLS <- c("12-3_P_DMSO_S1", "12-4_P_DMSO_S2")
  DIFF_COLS <- c("12-3_SS_DMSO_S5", "12-4_SS_DMSO_S6")
  BASE2 <- c("12-3_P_A196_S3", "12-4_P_A196_S4")
  DIFF2 <- c("12-3_SS_A196_S7", "12-4_SS_A196_S8")
} else if (ANALYSIS == "siCTCF") {
  # A196, SS
  COLS_MASK <- c("12-3_CI_SiControl_S21", "12-4_CI_SiControl_S22",
                 "12-3_CI_SiCTCF_S23", "12-4_CI_SiCTCF_S24")#,
  #"12-3_P_DMSO_S1", "12-4_P_DMSO_S2",
  #"12-3_P_A196_S3", "12-4_P_A196_S4")
  
  # special case where we want to sort by prolif.
  DIFF_COLS <- c("12-3_P_SiControl_S13", "12-4_P_SiControl_S14")
  BASE_COLS <- c("12-3_SS_SiControl_S17", "12-4_SS_SiControl_S18")
  DIFF2 <- c("12-3_P_SiCTCF_S15", "12-4_P_SiCTCF_S16")
  BASE2 <- c("12-3_SS_SiCTCF_S19", "12-4_SS_SiCTCF_S20")
} else {
  stop("Invalid analysis run")
}




if(GENE_SUBSET == "sig_in_q") {
  subset_deg_universe <- sigGenes.consensus
} else if (GENE_SUBSET == "down_in_q") {
  subset_deg_universe <- downGenes.consensus
} else if (GENE_SUBSET == "up_in_q") {
  subset_deg_universe <- upGenes.consensus
} else {
  stop("invalid OE run")
}
#-------------------------------------------------------------------------------
TRIPLICATE_DIR <- DESEQ_DIR
COUNTS.NORM.DEG.TRIPLICATE <-read.csv(TRIPLICATE_NORM_COUNTS_FILE) %>%
  dplyr::filter(gene_name %in% subset_deg_universe) # only the sig DEG
as_tibble(COUNTS.NORM.DEG.TRIPLICATE)


# remove gene names with multiple gene id matches
#   (also, the counts/normalization is probably messed up anyway)
COUNTS.NORM.DEG.TRIPLICATE.NONDUPLICATE_GENE_NAME <- COUNTS.NORM.DEG.TRIPLICATE %>%
  group_by(gene_name) %>%        # Group by gene_name
  dplyr::filter(dplyr::n() == 1) %>%           # Keep only groups with exactly one occurrence
  ungroup()
COUNTS.NORM.DEG.TRIPLICATE.NONDUPLICATE_GENE_NAME %>% as_tibble()

TRIPLICATE_IDS <- colnames(COUNTS.NORM.DEG.TRIPLICATE.NONDUPLICATE_GENE_NAME)[
  colnames(COUNTS.NORM.DEG.TRIPLICATE.NONDUPLICATE_GENE_NAME) != "gene_name"
]
TRIPLICATE_IDS

ref.RFP_OE <- file.path(DESEQ_RESULTS_DIR, "RFP_OE.sig_and_lfc",
                        "RFP_OE.sig_and_lfc.SS_vs_P.DEgenes.results.csv")
ref.DMSO   <- file.path(DESEQ_RESULTS_DIR, "DMSO.sig_and_lfc",
                        "DMSO.sig_and_lfc.SS_vs_P.DEgenes.results.csv")
ref.siCTCF <- file.path(DESEQ_RESULTS_DIR, "SiControl.sig_and_lfc",
                        "SiControl.sig_and_lfc.SS_vs_P.DEgenes.results.csv")

if (ANALYSIS == "OE") {
  REF.FILE <- ref.RFP_OE
} else if (ANALYSIS == "A196") {
  REF.FILE <- ref.DMSO
} else if (ANALYSIS == "siCTCF") {
  REF.FILE <- ref.siCTCF
}

deg.reference <- REF.FILE %>%
  read_csv() 
colnames(deg.reference)[1] <- 'gene_id'
deg.reference <- deg.reference %>%
  filter(padj < THRESHOLD_P) %>%
  filter(abs(log2FoldChange) > THRESHOLD_LFC) %>%
  right_join(geneID_to_geneName, by='gene_id') %>%
  group_by(gene_name) %>%        # Group by gene_name
  dplyr::filter(dplyr::n() == 1) %>%           # Keep only groups with exactly one occurrence
  ungroup() %>%
  dplyr::select(-gene_id)
deg.reference

if(GENE_SUBSET == "sig_in_q") {
  deg.reference <- deg.reference %>%
    arrange(log2FoldChange)
} else if (GENE_SUBSET == "down_in_q") {
  deg.reference <- deg.reference %>% filter(log2FoldChange < 0) %>%
    arrange(log2FoldChange)
} else if (GENE_SUBSET == "up_in_q") {
  deg.reference <- deg.reference %>% filter(log2FoldChange > 0) %>%
    arrange(desc(log2FoldChange))
} else {
  stop("invalid OE run")
}
deg.reference

if (GENE_UNIVERSE == 'deg_reference') {
  REF_TABLE <- deg.reference %>% dplyr::select(gene_name)
} else if (GENE_UNIVERSE == 'triplicate') {
  REF_TABLE <- COUNTS.NORM.DEG.TRIPLICATE.NONDUPLICATE_GENE_NAME
}

#-------------------------------------------------------------------------------
#COLS_MASK <- c()
#"CI_10.6_quant", "CI_12.3_quant", "CI_12.4_quant")
COLS_SELECT <- rownames(sampleTable)[!(rownames(sampleTable) %in% COLS_MASK)]
COLS_SELECT


COLS_QUIESCENT <- c(TRIPLICATE_IDS, COLS_SELECT)[grep("(^SS_)|(_SS_)|(^CI_)|(_CI_)", c(TRIPLICATE_IDS, COLS_SELECT))]
COLS_QUIESCENT

filter_na <- function(counts.mat) {
  rowSums(is.na(counts.mat)) == 0
}
clip <- function(m) {
  m[m < -3] <- -3
  m[m > 3] <- 3
  m
}
cal_z_score <- function(x, clipped=F){
  z <- t(scale(t(x))) #(x - mean(x)) / sd(x)
  
  if (!clipped) {
    return(z)
  }
  return(clip(z))
}
# sorting rules for the counts matrix (better heatmap visualization than clustering the rows)
sort_counts.by_refList <- function(counts_matrix, col_by="gene_name", refList=deg.reference$gene_name) {
  refList.df <- data.frame(refList)
  colnames(refList.df) <- col_by
  filtered_matrix <- counts_matrix[rownames(counts_matrix) %in% refList, , drop = FALSE]
  ordered_matrix <- filtered_matrix[match(refList, rownames(filtered_matrix)), , drop = FALSE]
  ordered_matrix <- ordered_matrix[rowSums(is.na(ordered_matrix)) == 0, , drop = FALSE]
  return(
    filtered_matrix
  )
}
sort_counts.by_zscore <- function(counts_matrix, 
                                     calc_cols=COLS_SELECT,
                                     base_cols=BASE_COLS,
                                     diff_cols=DIFF_COLS#, 
                                     #c("12-3_SS_CTCF_OE_S31", "12-4_SS_CTCF_OE_S32")
                                     #"12-3_CI_CTCF_OE_S35", "12-4_CI_CTCF_OE_S36")
                                     #diff_cols=c("12-3_CI_RFP_OE_S33", "12-4_CI_RFP_OE_S34", "12-3_SS_RFP_OE_S29", "12-4_SS_RFP_OE_S30")
) {
  # Ensure the input is a matrix
  if (!is.matrix(counts_matrix)) {
    stop("Input must be a matrix.")
  }
  
  # Calculate the sorting criterion for each row (z-score --> diff in rowMeans of blocked conditions)
  counts_matrix_zscore <- cal_z_score(counts_matrix[,calc_cols], clipped=T)
  # sorting_criterion <- rowMeans(counts_matrix_zscore[, diff_cols]) - rowMeans(counts_matrix_zscore[, base_cols]) /
  # abs(rowMeans(counts_matrix_zscore[, base_cols]))
  sorting_criterion <- rowMeans(counts_matrix_zscore[, diff_cols])
  
  # Sort the matrix based on the criterion
  if ( (!(ANALYSIS %in% c('OE', 'A196')) && GENE_SUBSET == "up_in_q") || ((ANALYSIS %in% c('OE', 'A196')) && GENE_SUBSET == "down_in_q") ) {
    counts_matrix <- counts_matrix[rev(order(sorting_criterion)), ]
  } else {
    counts_matrix <- counts_matrix[order(sorting_criterion), ] # rev()
  }
  return(counts_matrix)
}

#-------------------------------------------------------------------------------

# remove gene names with multiple gene id matches
#   (also, the counts/normalization is probably messed up anyway)
counts.norm.comparison <- counts(dds, normalized=T) %>%
  #.[, COLS_SELECT] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var='gene_id') %>%
  right_join(geneID_to_geneName, by='gene_id')
as_tibble(counts.norm.comparison)

counts.norm.comparison.duplicates <- counts.norm.comparison %>%
  group_by(gene_name) %>%        # Group by gene_name
  dplyr::filter(dplyr::n() > 1) %>%           # Keep only groups with exactly one occurrence
  ungroup() %>%
  dplyr::select(-gene_id)
counts.norm.comparison.duplicates$gene_name # make sure these aren't important. most should be RF[0-9]+ genes

counts.norm.comparison.nonduplicate_gene_name <- counts.norm.comparison %>%
  group_by(gene_name) %>%        # Group by gene_name
  dplyr::filter(dplyr::n() == 1) %>%           # Keep only groups with exactly one occurrence
  ungroup() %>%
  dplyr::select(-gene_id)
as_tibble(counts.norm.comparison.nonduplicate_gene_name)

# merge triplicate data with OE data
counts.norm.comparison.merged.df <- dplyr::left_join(
  REF_TABLE, #COUNTS.NORM.DEG.TRIPLICATE.NONDUPLICATE_GENE_NAME, 
  counts.norm.comparison.nonduplicate_gene_name,
  by="gene_name") %>%
  .[filter_na(.),]
counts.norm.comparison.merged <- counts.norm.comparison.merged.df %>%
  column_to_rownames(var="gene_name") %>%
  as.matrix()

counts.ref.all <- COUNTS.NORM.DEG.TRIPLICATE.NONDUPLICATE_GENE_NAME %>%
  column_to_rownames(var="gene_name") %>%
  as.matrix()

gene_intersect <- intersect(
  rownames(counts.norm.comparison.merged),
  COUNTS.NORM.DEG.TRIPLICATE.NONDUPLICATE_GENE_NAME$gene_name
)

hits <- length(unique(gene_intersect))
tot_available <- length(
  unique(rownames(counts.norm.comparison.merged))
)
tot_triplicate <- length(unique(
  COUNTS.NORM.DEG.TRIPLICATE.NONDUPLICATE_GENE_NAME$gene_name
))

hits
hits/tot_available
hits/tot_triplicate



if (GENE_UNIVERSE == "deg_reference") {
  counts.norm.comparison.merged <- counts.norm.comparison.merged %>%
    #sort_counts.by_refList()
    sort_counts.by_refList(refList = gene_intersect)
}


counts.norm.comparison.merged.sorted <- counts.norm.comparison.merged %>%
  sort_counts.by_zscore(base_cols=BASE2, diff_cols=DIFF2)



make_heatmap.comparison <- function(counts.mat=counts.norm.comparison.merged, #counts(dds, normalized=T),
                            show_rows=TRUE, cluster_columns=FALSE, cluster_rows=TRUE,
                            show_row_names=F,
                            filter_method=filter_na) {
  #COLORS <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  # COLORS <- colorRamp2(c(-3, -2, -1, 0, 1, 2, 3), c(COLORS[1], COLORS[10], COLORS[25], COLORS[45], COLORS[75], COLORS[90], COLORS[100]))
  # colorRamp2(c(-3, 0, 3), c("darkcyan", "white", "darkred")) #c("#1984c5", "white", "#e14b31"))
  COLORS <- colorRamp2(c(-3, -2, -1, 0, 1, 2, 3), c("#22CCFF", "#1A99C0", "#0B4455", "#000000", "#450005", "#9C000C", "#D0000F")) #c("#22CCFF", "#1A99C0", "#1C5566", "#000000", "#561116", "#9C000C", "#D0000F")) #c("#0099FF", "#0088CC", "#003355", "#000000", "#450005", "#AD000D", "#D0000F")) #c("#0099FF", "#000000", "#D0000F")) # #0066CC
  COLORS <- colorRamp2(c(-3, -2, -1, 0, 1, 2, 3), c("#22CCFF", "#1A99C0", "#0B4455", "#000000", "#450005", "#9C000C", "#AC000C"))
  l = filter_method(counts.mat)
  l[is.na(l)] = FALSE
  if(sum(l) == 0) return(NULL)
  
  m = counts.mat #counts(dds, normalized = TRUE)
  m = m[l, ]
  
  env$row_index = which(l)
  
  ht = Heatmap(cal_z_score(m, clipped=T), name = "z-score", heatmap_legend_param = list(at = c(-3, -2, -1, 0, 1, 2, 3)), col = COLORS,
               # top_annotation = HeatmapAnnotation(
               #   dex = colData(dds)$dex,
               #   sizeFactor = anno_points( estimateSizeFactorsForMatrix(counts(dds)) ) #anno_points(colData(dds)$sizeFactor)
               # ),
               cluster_columns = cluster_columns, cluster_rows = cluster_rows,
               show_row_names = show_row_names, show_column_names = TRUE, #row_km = 2,
               column_title = paste0(sum(l), " genes"),
               show_row_dend = show_rows, 
               use_raster=FALSE,
               column_dend_height = unit(50, "mm"),
               row_dend_width = unit(30, "mm")
  ) #+ 
  # Heatmap(log10(res$baseMean[l]+1), show_row_names = FALSE, width = unit(5, "mm"),
  #         name = "log10(baseMean+1)", show_column_names = FALSE) +
  # Heatmap(res$log2FoldChange[l], show_row_names = FALSE, width = unit(5, "mm"),
  #         name = "log2FoldChange", show_column_names = FALSE,
  #         col = colorRamp2(c(-2, 0, 2), c("green", "black", "red")))
  ht = draw(ht, merge_legend = TRUE)
  ht
}

ht.ref.all <- make_heatmap.comparison(
  counts.mat=counts.ref.all
)

if (GENE_UNIVERSE == 'triplicate') {
  ht.ref.clustered <- make_heatmap.comparison(
    counts.mat=counts.norm.comparison.merged.sorted[, TRIPLICATE_IDS], 
    cluster_rows=T, show_rows=F
  )
  ht.ref <- make_heatmap.comparison(
    counts.mat=counts.norm.comparison.merged[, TRIPLICATE_IDS], 
    cluster_rows=F, show_rows=F
  )
  
  ht.ref.sorted <- make_heatmap.comparison(
    counts.mat=counts.norm.comparison.merged.sorted[, TRIPLICATE_IDS], 
    cluster_rows=F, show_rows=F
  )
}
ht.all <- make_heatmap.comparison(counts.mat=counts.norm.comparison.merged.sorted, #filter_method=filter_OE, 
                          cluster_rows=F, show_rows=F)

ht.cols_select <- make_heatmap.comparison(
  counts.mat=counts.norm.comparison.merged[, COLS_SELECT],
  cluster_rows=F, show_rows=F
)

ht.cols_select.sorted <- make_heatmap.comparison(
  counts.mat=counts.norm.comparison.merged.sorted[, COLS_SELECT],
  #filter_method=filter_OE,
  cluster_rows=F, show_rows=F
)

f_out <- file.path(ANALYSIS_DIR, "heatmap_figures",
                   paste0(ANALYSIS, "_", GENE_SUBSET, ".tiff"))

if (GENE_UNIVERSE == 'deg_reference') {
  tiff(f_out, height=6890, width=5450, res=1200)
  print(ht.cols_select.sorted)
  dev.off()
}
