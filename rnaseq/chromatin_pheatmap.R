# Chromatin remodeler pheatmap — Figure 1B.
# This script is largely self-contained: it reads normalized count CSVs and
# the consensus DEG list from disk. Run after run_deseq2.R has produced
# outputs, or point config.R paths at pre-existing result files.
# Requires config.R variables.

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
#-------------------------------------------------------------------------------
# visualize direction of change for all significant DEGs from this Gene Ontology category: |log fold change| > 0. 
# Can increase to |log fold change| > 0.5 or 1 if narrowing down candidate genes
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD <- 0 

SIG_FILTER <- "all" # "all" differentially expressed genes, "up" (upregulated), or "down" (downregulated)

CLUSTER_ROWS <- ifelse(SIG_FILTER == "all", F, F)
#-------------------------------------------------------------------------------

# load files -------------------------------------------------------------------

# chromatin genes
chromatin_remodelers.file <- CHROMATIN_REMODELERS_FILE
setwd(DESEQ_DIR)
counts <- read.csv(NORM_COUNTS_FILE)

chromatin_remodelers <- read.table(chromatin_remodelers.file) %>% pull(V1) %>% sort()
chromatin_remodelers %>% head()
cohesins <- c("RAD21", "REC8", "SMC1A", "SMC3")
condensins <- c("SMC2", "SMC4", "NCAPD2", "NCAPD3", "NCAPG", "NCAPG2",
                "NCAPH", "NCAPH2")



gene_list <- list(`cohesin`=cohesins, `condensin`=condensins,
                  `chromatin remodeling`=chromatin_remodelers)
lengths(gene_list)


filter_deg <- function(deg_df, filter_name=SIG_FILTER) {
  # # need to fix the sig filter (abs() is bad (erroneous) approach)
  if (filter_name == "all") {
    return(
      deg_df %>%
        filter(
            (SS_P_log2FoldChange > LFC_THRESHOLD & CI_P_log2FoldChange > LFC_THRESHOLD) | 
              (SS_P_log2FoldChange < -LFC_THRESHOLD & CI_P_log2FoldChange < -LFC_THRESHOLD) &
              (SS_P_padj < PADJ_THRESHOLD & CI_P_padj < PADJ_THRESHOLD)
            )
        )
  }
  if (filter_name == "up") {
    return(
      deg_df %>%
        filter((SS_P_log2FoldChange > LFC_THRESHOLD & CI_P_log2FoldChange > LFC_THRESHOLD)) %>%
          filter(SS_P_padj < PADJ_THRESHOLD & CI_P_padj < PADJ_THRESHOLD)
    )
  }
  if (filter_name == "down") {
    return(
      deg_df %>%
        filter((SS_P_log2FoldChange < -LFC_THRESHOLD & CI_P_log2FoldChange < -LFC_THRESHOLD)) %>%
          filter(SS_P_padj < PADJ_THRESHOLD & CI_P_padj < PADJ_THRESHOLD)
    )
  }
  stop("undefined filter name")
}

deseq.sigGenes.df <- CONSENSUS_DEG_FILE %>%
  read_csv() %>%
  filter_deg()# already DEGs, but sometimes we want additional stringent LFC threshold filter, etc.
deseq.sigGenes <- deseq.sigGenes.df %>%
  pull(gene_name)
deseq.sigGenes %>% head()
chromatin_intersect <- deseq.sigGenes.df %>%
  filter( 
    (gene_name %in% unlist(gene_list)) 
  )
dim(deseq.sigGenes.df)
dim(chromatin_intersect)
dim(chromatin_intersect %>% filter(SS_P_log2FoldChange > 0))
dim(chromatin_intersect %>% filter(SS_P_log2FoldChange < 0))


deseq.allGenes.df <- read_csv(NORM_COUNTS_FILE)

length(chromatin_remodelers)
length(intersect(chromatin_remodelers, deseq.allGenes.df$gene_name))

# filtering counts file ---------------------------------------------------
counts_filter <- counts %>% filter( 
  (gene_name %in% unlist(gene_list) & gene_name %in% deseq.sigGenes)
) %>%
  distinct()

gene_list.subset <- list()
gene_subset <- c()
for (n in names(gene_list)) {
  v <- gene_list[[n]]
  s <- v[v %in% counts_filter$gene_name]
  gene_list.subset[[n]] <- s
  gene_subset <- c(gene_subset, s)
}

gene_order <- data.frame(gene_name = gene_subset, order = seq_along(gene_subset))
# Reorder the matrix using left_join
counts_filter <- gene_order %>%
  left_join(counts_filter, by = "gene_name") %>% # Join to bring matrix rows into order
  arrange(order) %>%                        # Arrange rows by the order column
  dplyr::select(-order) 

stopifnot(all(counts_filter$gene_name == gene_subset) & all(counts_filter$gene_name %in% unlist(gene_list)))
head(counts_filter)
counts_filter_rownames<-counts_filter %>% remove_rownames %>% column_to_rownames(var="gene_name")
head(counts_filter_rownames)

# row colors and groups
annotation_row <- data.frame(`gene`=factor(rep(names(gene_list.subset), lengths(gene_list.subset))))
row.names(annotation_row) <- row.names(counts_filter_rownames)
annotation_row

# color scheme of the heatmap --------------------------
n_colors <- 100
color_range = c("#1A99C0", "#000000", "#AC000C") #create a color range
color_range_expand = colorRampPalette(color_range)(n_colors) #expansion of the color range

# calculate z score  -----------------------------------
display.brewer.all()
color_brewer = brewer.pal(n=11, name="RdBu")
color_brewer_expand = colorRampPalette(color_brewer)(50)
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
anchor_values <- function(mat, min_val, max_val){
  for (r in 1:(length(mat[,1])-1)) {
    for (c in 1:(length(mat[r,])-1)) {
      mat[r, c] <- min(max_val, mat[r, c]) %>% max(min_val)
    }
  }
  mat
}
counts_filter_rownames_zscore <- t(apply(counts_filter_rownames, 1, cal_z_score))

# add column bars for the conditions --------------------------------------
sample_col<-read.csv("sample_col.csv")
row.names(sample_col)<-colnames(counts_filter_rownames_zscore)
head(sample_col)


# pheatmap --------------------------------------
library(pheatmap)
library(grid)
library(gtable)
library(ggplot2)

plot.out <- pheatmap(
  counts_filter_rownames_zscore %>% anchor_values(-3, 3), 
  color = color_range_expand,
  cellheight = 10, cellwidth = 50,
  cluster_cols = FALSE, cluster_rows = CLUSTER_ROWS, 
  annotation_col = sample_col,
  annotation_row = annotation_row,
  annotation_colors = list(
    sample = c(P = 'coral1', CI = 'pink', SS = 'dodgerblue', CIR = 'yellowgreen', SSR = 'darkgoldenrod'),
    gene = c(cohesin = 'skyblue', condensin = 'darkgreen', `chromatin remodeling` = 'darkgray')
  ),
  breaks = seq(-3, 3, length.out = n_colors + 1),
  legend_breaks = c(-3, -2, -1, 0, 1, 2, 3),
  annotation_names_row = TRUE
)

tiff(
  filename = paste0("deseq_chromatin.lfc", LFC_THRESHOLD, ".", SIG_FILTER, "final.tiff"),
  width = 40, height = 40, units = "cm", res = 1200
)
print(plot.out)
dev.off()


anno_df <- annotation_row %>%
  mutate(
    gene_label = rownames(annotation_row),  # keep rownames in a new column
    .order = match(rownames(annotation_row), rownames(counts_filter_rownames_zscore))
  )
anno_plot <- ggplot(anno_df,
                    aes(x = 1, y = reorder(gene_label, -.order), fill = gene)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(
    cohesin = "skyblue",
    condensin = "darkgreen",
    `chromatin remodeling` = "darkgray"
  )) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 2, 0, 0),
    axis.text.y = element_text(size = 7, hjust = 1)
  )
png(
  filename = paste0("deseq_chromatin.lfc", LFC_THRESHOLD, ".", SIG_FILTER, "final.anno_row.png"),
  width = 40, height = 40, units = "cm", res = 1200
)
print(anno_plot)
dev.off()

counts_filter_rownames_zscore[, 2] %>% hist()

# inspect the pheatmap gene clusters (reruns clustering step) ------------------
hc <- hclust(dist(counts_filter_rownames_zscore))
plot(hc, cex=0.6)
k <- 2
hc_groups <- cutree(hc, k=k)
hc_groups %>% head()

analysis_name <- paste0("deseq_chromatin.lfc", LFC_THRESHOLD,".", SIG_FILTER)
genes_by_group <- list()
for (i in 1:k) {
  genes_by_group[[i]] <- names(hc_groups)[hc_groups == i]
  if (SIG_FILTER == 'all') {
    writeLines(genes_by_group[[i]], paste0(analysis_name, ".chromatin.group.", i, ".txt"))
    
    counts_by_group <- counts %>%
      filter(gene_name %in% genes_by_group[[i]])
    write_csv(counts_by_group, paste0(analysis_name, ".chromatin.group.", i, ".norm_log2_counts_final_noquant.csv"))
  }
}
genes_by_group

counts_chromatin <- counts %>% filter(gene_name %in% chromatin_intersect$gene_name)
counts_chromatin
write_csv(counts_chromatin, paste0(analysis_name, ".chromatin", ".norm_log2_counts_final_noquant.csv"))