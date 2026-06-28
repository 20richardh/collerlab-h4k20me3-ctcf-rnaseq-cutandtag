# Core DESeq2 / tximport workflow.
# Sourced by run_deseq2.R after config.R and a design file have been loaded.
# Do not run this file directly.

library(tidyverse)

DATA_COLS <- c("bioReplicate", "quiescence", "cell_type", "gene_expression", "is_OE")

#creating the tx2gene data frame containing transcript id and gene id
library(AnnotationDbi)
library(GenomicFeatures)
txdb<- txdbmaker::makeTxDbFromGFF(file=genome_path)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)

#counting the gene counts from transcripts abundance (with no offset)
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)
names(txi)
head(txi$counts)

#for no offset counts, DESeqDataSetFromTximport should be used first
#sampleTable <- data.frame(condition = factor(rep(c("P", "CI"), each = 3)))
sampleTable <- data.frame(condition = samples[,compare_by$col[1]], row.names=samples$run)
#rownames(sampleTable) <- colnames(txi$counts)
sampleTable
stopifnot(all(colnames(txi$counts) == rownames(sampleTable)))
print("NOTE: Manually check that the condition matches the row name (samples$run) in the sampleTable")

library(DESeq2)
dds.pre <- DESeqDataSetFromTximport(txi, sampleTable, ~condition) #used this for DESeq analysis below
dds.pre$condition <- relevel(dds.pre$condition, ref = compare_by$ref[1])#"P") #for log2(quiescence protocol/P)

####DESeq2 analysis using the dds.pre obtained from DESeqDataSetFromTximport
# keep <- rowSums(counts(dds.pre)) >= 10 #sum of total counts >=10
print("NOTE: Check that the smallest group size cutoff below is accurate")
group_sizes <- dplyr::count(sampleTable, condition)
group_sizes.smallest <- min(group_sizes$n)
group_sizes.smallest
keep <- rowSums(counts(dds.pre) >= 10) >= group_sizes.smallest
dds.pre <- dds.pre[keep,]
dds <- DESeq(dds.pre)
res <- results(dds) #this will give log2(FC)
head(res)
summary(res)

# summary_filepath <- paste(analysis_name, 'results.summary.txt', sep='.')
# fileConn <- file(summary_filepath)
# string_data <- summary(res) %>% toString()
# writeLines(string_data, fileConn)
# close(fileConn)

# INSPECT SOME INTERESTING OR SUSPECT GENES. ALSO VISUALIZE WITH DESEQ2 PLOTCOUNTS
genes.lowcounts <- rownames(res[is.na(res$padj),])
head(genes.lowcounts, n=25)
counts.lowcounts <- counts(dds.pre)[rownames(counts(dds.pre)) %in% genes.lowcounts,]
#head(counts.lowcounts)

if (length(genes.lowcounts)) {
  target_gene <- genes.lowcounts[[1]]
  plotCounts(dds.pre, target_gene)
} else {
  print("No low-count genes based on DESeq2 dispersion")
}

# counts
par(mar=c(8,5,2,2))
plot.cooks <- boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
plot.cooks
tiff(paste(analysis_name, 'cook.tiff', sep='.'))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

# dispersion. black is high dispersion, blue circle is within the curve (supposedly false positive) 
#   (https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html)
plotDispEsts(dds)
tiff(paste(analysis_name, 'dispersion.tiff', sep='.'))
plotDispEsts(dds)
dev.off()



##Log fold change shrinkage for visualization and ranking
library(apeglm)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef=paste("condition", comparison_name, sep='_'), type="ashr")
head(resLFC) #this will give log2(P/CI)

##MA-plot
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))


##print files
# Save tximport counts (raw counts from salmon quantification)
results_file <- paste(analysis_name, "tximport_counts.csv", sep='.')
write.csv(txi$counts, file = results_file, row.names = TRUE)

# Get log2 normalized counts from DESeq2
log2_normalized_counts <- log2(counts(dds, normalized = TRUE) + 1)

# Save log2 normalized counts
results_file <- paste(analysis_name, "log2_normalized_counts.csv", sep='.')
write.csv(log2_normalized_counts, file = results_file, row.names = TRUE)

results_file <- paste(analysis_name, "DEgenes.results.csv", sep='.')
results_file.LFC <- paste(analysis_name, "DEgenes.results.lfcshrink.csv", sep='.')
write.csv(as.data.frame(res), 
          file=results_file)

write.csv(as.data.frame(resLFC), 
          file=results_file.LFC)

##data transformation by variance stabilizing transformations (VST) method
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)




#####To add the gene names as obtained from gencode v29 annotation (all genes)
#Add "gene_id" to first column of the deseq2 output file

##load packages
#install.packages("tidyverse")
#library(tidyverse)

##load tables for gencode gtf (converted to csv) and DESeq2 results
gencode<-read_csv(genome_csv_path)  %>% distinct()
deseq<-read_csv(results_file)
colnames(deseq)[1] <- 'gene_id'
as_tibble(deseq)

##prepare the gencode csv file for further analysis
gencode_select <-gencode # %>% select(gene_id, gene_name, gene_type)
gencode_unique<-distinct(gencode_select)
geneID_to_geneName <- dplyr::select(gencode, gene_id, gene_name) %>% distinct()

##add the gene names to DESeq2 results
deseq_with_gene_name<-inner_join(deseq, gencode_unique, by="gene_id") 

##select significant genes only (unstranded)
#deseq_lincRNA <- deseq_with_gene_name %>% dplyr::filter(grepl("lincRNA", gene_type)) 
deseq_genes_sig <- deseq_with_gene_name %>% dplyr::filter(padj<THRESHOLD_P) 

##write files
allgenes_file <- paste(analysis_name, "gencodev29.allgenes.csv", sep='.')
allgenes_sig_file <- paste(analysis_name, "gencodev29.allgenes.sig.csv", sep='.')
write_csv(deseq_with_gene_name, allgenes_file)
write_csv(deseq_genes_sig, allgenes_sig_file)




library(ggfortify)
library(PCAtools)
library(viridis)
library(ggplot2)
library(grid)
library(gridExtra)

# ---------------------------------------------------------------------------
# Build a single biplot (REMOVED LEGEND, saved separately)
# ---------------------------------------------------------------------------
make_biplot <- function(pca_data, x, y,
                        colby, colkey,
                        shape_var, shapekey,
                        pointSize = 4,
                        xlims = c(-30, 30),
                        ylims = c(-30, 30)) {
  
  p <- biplot(pca_data,
              x = x, y = y,
              colby        = colby,     colkey    = colkey,
              shape        = shape_var, shapekey  = shapekey,
              pointSize    = pointSize,
              showLoadings = FALSE,
              lab          = NULL) +
    xlim(xlims) + ylim(ylims) +
    theme_classic() +
    theme(legend.position = "none")   # legend saved separately
  
  # Replace underscores with spaces in discrete scale labels
  p$scales$scales <- lapply(p$scales$scales, function(s) {
    if (inherits(s, "ScaleDiscrete") && !is.null(s$labels)) {
      s$labels <- gsub("_", " ", s$labels)
    }
    s
  })
  
  return(p)
}

# ---------------------------------------------------------------------------
# extract the legend grob from a ggplot object.
# ---------------------------------------------------------------------------
extract_legend <- function(p) {
  gt  <- ggplot_gtable(ggplot_build(p))
  # "guide-box" may appear as a single grob or as multiple — grab all matches
  idx <- grep("guide-box", gt$layout$name)
  if (length(idx) == 0) stop("No legend found in plot — check legend.position is not 'none'.")
  # If multiple, cowplot-style: wrap them together
  if (length(idx) == 1) return(gt$grobs[[idx]])
  do.call(gridExtra::arrangeGrob, gt$grobs[idx])
}

# ---------------------------------------------------------------------------
# Draw legend for the biplots that don't have replicates labeled
#   colour aesthetic → shown as a horizontal dash (linetype segment)
#   shape  aesthetic → shown as a point in black
# ---------------------------------------------------------------------------
make_legend_standard <- function(col_colors, shape_filled,
                                 color_by_label, shape_label = "Cell type") {
  df <- expand.grid(
    color_by  = factor(names(col_colors),  levels = names(col_colors)),
    cell_type = factor(names(shape_filled), levels = names(shape_filled)),
    stringsAsFactors = FALSE
  )
  df$x <- seq_len(nrow(df))
  df$y <- 0
  
  p_leg <- ggplot(df, aes(x = x, y = y,
                          colour = color_by,
                          shape  = cell_type)) +
    geom_point(size = 4) +
    scale_colour_manual(
      name   = gsub("_", " ", color_by_label),
      values = col_colors,
      labels = gsub("_", " ", names(col_colors))
    ) +
    scale_shape_manual(
      name   = gsub("_", " ", shape_label),
      values = shape_filled,
      labels = gsub("_", " ", names(shape_filled))
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(shape     = NA,
                            linetype  = "solid",
                            linewidth = 1.5,
                            size      = 6)
      ),
      shape = guide_legend(
        override.aes = list(colour = "black",
                            size   = 3)
      )
    ) +
    theme_void() +
    theme(legend.position = "right",
          legend.title    = element_text(size = 11, face = "bold"),
          legend.text     = element_text(size = 10))
  
  extract_legend(p_leg)
}

# ---------------------------------------------------------------------------
# Extract replicate ID from sample name (e.g. "12-3_xxx" -> "12-3")
# ---------------------------------------------------------------------------
extract_replicate <- function(sample_names) {
  sub("^([^_]+)_.*$", "\\1", sample_names)
}

# ---------------------------------------------------------------------------
# Main PCA function
# ---------------------------------------------------------------------------
deseq2_pca <- function(vsd, ntop = 2000, sample_metadata = samples,
                       returnData = FALSE,
                       color_by   = list(idx  = "gene_expression_idx",
                                         name = "gene_expression"),
                       df.info    = gencode_unique,
                       merge_by   = "gene_id") {
  
  sample_metadata           <- as.data.frame(sample_metadata)
  rownames(sample_metadata) <- sample_metadata$run
  # "color_by" column holds the actual values used for colouring
  sample_metadata$color_by  <- sample_metadata[[color_by$name]]
  
  # --- variance-based gene selection ---
  rv        <- rowVars(assay(vsd))
  select    <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  vsd_assay <- assay(vsd)[select, ]
  
  # --- PCAtools PCA ---
  pca_data <- pca(vsd_assay, metadata = sample_metadata)
  
  rotation.df               <- pca_data$loadings
  rotation.df[["gene_id"]]  <- rownames(rotation.df)
  
  # --- Aesthetics ---------------------------------------------------------
  color_levels <- sort(unique(sample_metadata$color_by))
  if (compare_by$perturb == "CTCF_OE") {
    col_colors   <- setNames(c('darkgreen', 'black'), color_levels)
  } else if (compare_by$perturb == "A196") {
    col_colors <- setNames(c('magenta', 'black'), color_levels)
  } else {
    col_colors <- setNames(c('black', 'red'), color_levels)
  }
  print(col_colors)
  # Filled symbols for rep 1; open equivalents for rep 2
  shape_filled <- c('P' = 16, 'SS' = 17)   # filled circle, filled triangle
  shape_open   <- c('P' =  1, 'SS' =  2)   # open  circle, open  triangle
  
  # --- Replicate metadata -------------------------------------------------
  pca_data$metadata$replicate <- extract_replicate(rownames(pca_data$metadata))
  rep_levels_sorted <- sort(unique(pca_data$metadata$replicate))  # e.g. c("12-3","12-4")
  
  # ---------------------------------------------------------------------------
  # Convenience wrappers
  # ---------------------------------------------------------------------------
  make_bp <- function(x, y, xlims, ylims)
    make_biplot(pca_data,
                x         = x,          y         = y,
                colby     = "color_by", colkey    = col_colors,
                shape_var = "cell_type", shapekey  = shape_filled,
                pointSize = 4,
                xlims     = xlims,      ylims     = ylims)
  
  make_bp_rep <- function(x, y, xlims, ylims)
    make_biplot_rep(pca_data,
                    x              = x,                    y              = y,
                    color_by_col   = "color_by",           col_colors     = col_colors,
                    shape_col      = "cell_type",
                    shape_filled   = shape_filled,         shape_open     = shape_open,
                    rep_col        = "replicate",          rep_levels     = rep_levels_sorted,
                    color_by_label = color_by$name,
                    xlims          = xlims,                ylims          = ylims)
  
  # ---------------------------------------------------------------------------
  # All axis-limit variants — every unique x/y combination from the original
  # script (active and commented-out lines), named x{xlim}y{ylim}.
  # ---------------------------------------------------------------------------
  
  # PC1 vs PC2:  x±30 y±30  |  x±30 y±50  |  x±50 y±50
  plot.1.2.x30y30     <- make_bp    ("PC1","PC2", c(-30,30), c(-30,30))
  plot.1.2.x30y50     <- make_bp    ("PC1","PC2", c(-30,30), c(-50,50))
  plot.1.2.x50y50     <- make_bp    ("PC1","PC2", c(-50,50), c(-50,50))
  plot.1.2.x30y30.rep <- make_bp_rep("PC1","PC2", c(-30,30), c(-30,30))
  plot.1.2.x30y50.rep <- make_bp_rep("PC1","PC2", c(-30,30), c(-50,50))
  plot.1.2.x50y50.rep <- make_bp_rep("PC1","PC2", c(-50,50), c(-50,50))
  
  # PC1 vs PC3:  x±30 y±30  |  x±50 y±50
  plot.1.3.x30y30     <- make_bp    ("PC1","PC3", c(-30,30), c(-30,30))
  plot.1.3.x50y50     <- make_bp    ("PC1","PC3", c(-50,50), c(-50,50))
  plot.1.3.x30y30.rep <- make_bp_rep("PC1","PC3", c(-30,30), c(-30,30))
  plot.1.3.x50y50.rep <- make_bp_rep("PC1","PC3", c(-50,50), c(-50,50))
  
  # PC2 vs PC3:  x±30 y±30  |  x±50 y±50
  plot.2.3.x30y30     <- make_bp    ("PC2","PC3", c(-30,30), c(-30,30))
  plot.2.3.x50y50     <- make_bp    ("PC2","PC3", c(-50,50), c(-50,50))
  plot.2.3.x30y30.rep <- make_bp_rep("PC2","PC3", c(-30,30), c(-30,30))
  plot.2.3.x50y50.rep <- make_bp_rep("PC2","PC3", c(-50,50), c(-50,50))
  
  # PC3 vs PC1:  x±20 y±35  |  x±30 y±35  |  x±50 y±50
  plot.3.1.x20y35     <- make_bp    ("PC3","PC1", c(-20,20), c(-35,35))
  plot.3.1.x30y35     <- make_bp    ("PC3","PC1", c(-30,30), c(-35,35))
  plot.3.1.x50y50     <- make_bp    ("PC3","PC1", c(-50,50), c(-50,50))
  plot.3.1.x20y35.rep <- make_bp_rep("PC3","PC1", c(-20,20), c(-35,35))
  plot.3.1.x30y35.rep <- make_bp_rep("PC3","PC1", c(-30,30), c(-35,35))
  plot.3.1.x50y50.rep <- make_bp_rep("PC3","PC1", c(-50,50), c(-50,50))
  
  # PC3 vs PC4:  x±30 y±30  |  x±50 y±50
  plot.3.4.x30y30     <- make_bp    ("PC3","PC4", c(-30,30), c(-30,30))
  plot.3.4.x50y50     <- make_bp    ("PC3","PC4", c(-50,50), c(-50,50))
  plot.3.4.x30y30.rep <- make_bp_rep("PC3","PC4", c(-30,30), c(-30,30))
  plot.3.4.x50y50.rep <- make_bp_rep("PC3","PC4", c(-50,50), c(-50,50))
  
  # Misc
  plot.screeplot <- screeplot(pca_data)
  plot.loadings  <- plotloadings(pca_data, labSize = 3) +
    theme(legend.position = "none")
  plot.pairsplot <- pairsplot(pca_data, colby = "color_by", colkey = col_colors)
  
  # --- Legends (built once per PCA run, saved separately) ----------------
  legend.standard <- make_legend_standard(
    col_colors     = col_colors,
    shape_filled   = shape_filled,
    color_by_label = color_by$name
  )
  
  # shape_map and shape_labels are needed for the rep legend;
  # reconstruct them here to keep make_biplot_rep self-contained
  shape_map_leg <- c(
    setNames(shape_filled, paste0(names(shape_filled), ".", rep_levels_sorted[1])),
    setNames(shape_open,   paste0(names(shape_open),   ".", rep_levels_sorted[2]))
  )
  shape_labels_leg <- setNames(
    paste0(sub("\\..*", "", names(shape_map_leg)),
           " (rep ", sub(".*\\.", "", names(shape_map_leg)), ")"),
    names(shape_map_leg)
  )
  legend.rep <- make_legend_rep(
    col_colors     = col_colors,
    shape_map      = shape_map_leg,
    shape_labels   = shape_labels_leg,
    color_by_label = color_by$name
  )
  
  # ---------------------------------------------------------------------------
  return(list(
    data = pca_data,
    # legends
    legend.standard     = legend.standard,
    legend.rep          = legend.rep,
    # PC1 vs PC2
    plot.1.2.x30y30     = plot.1.2.x30y30,
    plot.1.2.x30y50     = plot.1.2.x30y50,
    plot.1.2.x50y50     = plot.1.2.x50y50,
    plot.1.2.x30y30.rep = plot.1.2.x30y30.rep,
    plot.1.2.x30y50.rep = plot.1.2.x30y50.rep,
    plot.1.2.x50y50.rep = plot.1.2.x50y50.rep,
    # PC1 vs PC3
    plot.1.3.x30y30     = plot.1.3.x30y30,
    plot.1.3.x50y50     = plot.1.3.x50y50,
    plot.1.3.x30y30.rep = plot.1.3.x30y30.rep,
    plot.1.3.x50y50.rep = plot.1.3.x50y50.rep,
    # PC2 vs PC3
    plot.2.3.x30y30     = plot.2.3.x30y30,
    plot.2.3.x50y50     = plot.2.3.x50y50,
    plot.2.3.x30y30.rep = plot.2.3.x30y30.rep,
    plot.2.3.x50y50.rep = plot.2.3.x50y50.rep,
    # PC3 vs PC1
    plot.3.1.x20y35     = plot.3.1.x20y35,
    plot.3.1.x30y35     = plot.3.1.x30y35,
    plot.3.1.x50y50     = plot.3.1.x50y50,
    plot.3.1.x20y35.rep = plot.3.1.x20y35.rep,
    plot.3.1.x30y35.rep = plot.3.1.x30y35.rep,
    plot.3.1.x50y50.rep = plot.3.1.x50y50.rep,
    # PC3 vs PC4
    plot.3.4.x30y30     = plot.3.4.x30y30,
    plot.3.4.x50y50     = plot.3.4.x50y50,
    plot.3.4.x30y30.rep = plot.3.4.x30y30.rep,
    plot.3.4.x50y50.rep = plot.3.4.x50y50.rep,
    # misc
    plot.screeplot  = plot.screeplot,
    plot.loadings   = plot.loadings,
    plot.pairsplot  = plot.pairsplot,
    rotation        = left_join(rotation.df, df.info, by = merge_by)
  ))
}

# ---------------------------------------------------------------------------
# Save all plots
# ---------------------------------------------------------------------------
if (!SKIP_PCA) {
  for (ntop_val in c(500, 1000, 2000)) {
    pca.out <- deseq2_pca(vsd, ntop = ntop_val)
    
    # ntop embedded in every output filename, e.g. "analysis.pca.ntop500.PC1_PC2.lim30.tiff"
    pca_path <- function(...) paste(analysis_name, "pca", paste0("ntop", ntop_val), ..., sep = '.')
    
    # --- Legends (one per ntop run) ---
    tiff(pca_path("legend.standard.tiff"), width = 1200, height = 1200, res = 300)
    grid.draw(pca.out$legend.standard)
    dev.off()
    tiff(pca_path("legend.rep.tiff"), width = 1200, height = 1200, res = 300)
    grid.draw(pca.out$legend.rep)
    dev.off()
    
    # --- Misc ---
    tiff(pca_path("screeplot.tiff")); print(pca.out$plot.screeplot); dev.off()
    tiff(pca_path("loadings.tiff"));  print(pca.out$plot.loadings);  dev.off()
    tiff(pca_path("pairsplot.tiff")); print(pca.out$plot.pairsplot); dev.off()
    
    # --- All biplot variants ---
    biplot_variants <- list(
      # PC1 vs PC2
      list("PC1_PC2.x30y30",     pca.out$plot.1.2.x30y30),
      list("PC1_PC2.x30y50",     pca.out$plot.1.2.x30y50),
      list("PC1_PC2.x50y50",     pca.out$plot.1.2.x50y50),
      list("PC1_PC2.x30y30.rep", pca.out$plot.1.2.x30y30.rep),
      list("PC1_PC2.x30y50.rep", pca.out$plot.1.2.x30y50.rep),
      list("PC1_PC2.x50y50.rep", pca.out$plot.1.2.x50y50.rep),
      # PC1 vs PC3
      list("PC1_PC3.x30y30",     pca.out$plot.1.3.x30y30),
      list("PC1_PC3.x50y50",     pca.out$plot.1.3.x50y50),
      list("PC1_PC3.x30y30.rep", pca.out$plot.1.3.x30y30.rep),
      list("PC1_PC3.x50y50.rep", pca.out$plot.1.3.x50y50.rep),
      # PC2 vs PC3
      list("PC2_PC3.x30y30",     pca.out$plot.2.3.x30y30),
      list("PC2_PC3.x50y50",     pca.out$plot.2.3.x50y50),
      list("PC2_PC3.x30y30.rep", pca.out$plot.2.3.x30y30.rep),
      list("PC2_PC3.x50y50.rep", pca.out$plot.2.3.x50y50.rep),
      # PC3 vs PC1
      list("PC3_PC1.x20y35",     pca.out$plot.3.1.x20y35),
      list("PC3_PC1.x30y35",     pca.out$plot.3.1.x30y35),
      list("PC3_PC1.x50y50",     pca.out$plot.3.1.x50y50),
      list("PC3_PC1.x20y35.rep", pca.out$plot.3.1.x20y35.rep),
      list("PC3_PC1.x30y35.rep", pca.out$plot.3.1.x30y35.rep),
      list("PC3_PC1.x50y50.rep", pca.out$plot.3.1.x50y50.rep),
      # PC3 vs PC4
      list("PC3_PC4.x30y30",     pca.out$plot.3.4.x30y30),
      list("PC3_PC4.x50y50",     pca.out$plot.3.4.x50y50),
      list("PC3_PC4.x30y30.rep", pca.out$plot.3.4.x30y30.rep),
      list("PC3_PC4.x50y50.rep", pca.out$plot.3.4.x50y50.rep)
    )
    
    for (v in biplot_variants) {
      tiff(pca_path(v[[1]], "tiff"), width = 6000, height = 6000, res = 1200)
      print(v[[2]])
      dev.off()
    }
    
    # --- Rotation ---
    write_csv(pca.out$rotation, pca_path("rotation.csv"))
  }
}


# POISSON DISTANCE, WITH PHEATMAP. Which samples are most similar to each other?
library(PoiClaClu)
library(RColorBrewer)
library(pheatmap)

poisd.counts <- counts(dds.pre)
#poisd.counts <- counts(dds.pre)[rownames(counts(dds.pre)) %in% pca.top_rotation_genes.df$gene_id, ]
poisd <- PoissonDistance(t(poisd.counts))
samplePoisDistMatrix <- as.matrix( poisd$dd )
colnames(samplePoisDistMatrix) <- colnames(counts(dds.pre)) #NULL
rownames(samplePoisDistMatrix) <- colnames(counts(dds.pre))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

cell_dim <- 5*round(1/5*sqrt(10000/length(colnames(poisd.counts))))
cell_dim
poisson_heatmap <- pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd, clustering_distance_cols = poisd$dd, col = colors,
         main = "Poisson distances",
         cellheight = cell_dim, cellwidth=cell_dim)
poisson_heatmap

poisson_distance.tiff <- paste(analysis_name, "poisson_distance", "tiff", sep='.')
tiff(poisson_distance.tiff, width=1000, height=1000)
print(poisson_heatmap)
dev.off()

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
dist_heatmap <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dist_heatmap
distance.tiff <- paste(analysis_name, "distance", "tiff", sep='.')
tiff(distance.tiff, width=1000, height=1000)
print(dist_heatmap)
dev.off()


# Transform count data using the variance stablilizing transform
library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)
vsd_df <- assay(vsd) %>% as.data.frame()
vsd_df$gene_id <- rownames(vsd_df)
vsd_df <- left_join(vsd_df, gencode_unique, by='gene_id')
head(vsd_df)

# Keep only the significantly differentiated genes passing log2FoldChange threshold
#sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= THRESHOLD_P & abs(deseq2ResDF$log2FoldChange) > THRESHOLD_LFC,])
sigGenes <- deseq_genes_sig$gene_id
vsd_sig <- vsd_df[vsd_df$gene_id %in% sigGenes,]
head(vsd_sig)

get_top_genes <- function(deseq_genes_sig, n, upregulated) {
  # TODO: sort the genes (alphabetically?) OR cluster them (see OneNote on visualization)
  if (upregulated) {
    deseq_genes_sig <- deseq_genes_sig %>% filter(log2FoldChange > 0)
  } else {
    deseq_genes_sig <- deseq_genes_sig %>% filter(log2FoldChange < 0)
  }
  top_upregulated <- deseq_genes_sig$gene_id[head(order(deseq_genes_sig$log2FoldChange, decreasing=upregulated), n=n)]
  vsd_sig.top_upregulated <- vsd_sig[vsd_sig$gene_id %in% top_upregulated,]
  return(vsd_sig.top_upregulated)
}

n <- 50
#UPREGULATED
vsd_sig.top_upregulated <- get_top_genes(deseq_genes_sig, n, upregulated=TRUE)
as_tibble(vsd_sig.top_upregulated)
#DOWNREGULATED
vsd_sig.top_downregulated <- get_top_genes(deseq_genes_sig, n, upregulated=FALSE)
as_tibble(vsd_sig.top_downregulated)

# sanity check:
check_upregulated <- deseq_genes_sig[deseq_genes_sig$gene_id %in% vsd_sig.top_upregulated$gene_id, ]
as_tibble(check_upregulated)
stopifnot(length(check_upregulated[[1]]) <= n)
check_downregulated <- deseq_genes_sig[deseq_genes_sig$gene_id %in% vsd_sig.top_downregulated$gene_id, ]
as_tibble(check_downregulated)
stopifnot(length(check_downregulated[[1]]) <= n)


library(reshape2)


# ComplexHeatMap
#BiocManager::install("ComplexHeatmap")
#BiocManager::install("InteractiveComplexHeatmap")
library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
library(circlize)
library(GetoptLong)

env = new.env()
filter_by_sig <- function(res, base_mean = 0, log2fc = 0) {
  res$padj <= THRESHOLD_P & res$baseMean >= base_mean & 
    abs(res$log2FoldChange) >= log2fc
}
filter_by_lfc <- function(res, base_mean = 0, log2fc = 1) {
  res$padj <= THRESHOLD_P & res$baseMean >= base_mean & 
    abs(res$log2FoldChange) >= log2fc
}
make_heatmap <- function(counts.df=counts(dds, normalized=T),
                         show_rows=TRUE, cluster_columns=FALSE, cluster_rows=TRUE,
                         row_km=2,
                         filter_method=filter_by_sig) {
  l = filter_method(res)
  l[is.na(l)] = FALSE
    
    if(sum(l) == 0) return(NULL)
    
    m = counts.df #counts(dds, normalized = TRUE)
    m = m[l, ]
    
    env$row_index = which(l)
    
    ht = Heatmap(t(scale(t(m))), name = "z-score", col = colorRamp2(c(-3, 0, 3), c("#1984c5", "white", "#e14b31")),
                 top_annotation = HeatmapAnnotation(
                   dex = colData(dds)$dex,
                   sizeFactor = anno_points( estimateSizeFactorsForMatrix(counts(dds)) ) #anno_points(colData(dds)$sizeFactor)
                 ),
                 cluster_columns = cluster_columns, cluster_rows = cluster_rows,
                 show_row_names = FALSE, show_column_names = TRUE, row_km = row_km,
                 column_title = paste0(sum(l), " genes"),
                 show_row_dend = show_rows, 
                 use_raster=FALSE,
                 column_dend_height = unit(50, "mm"),
                 row_dend_width = unit(30, "mm")
                 ) + 
      Heatmap(log10(res$baseMean[l]+1), show_row_names = FALSE, width = unit(5, "mm"),
              name = "log10(baseMean+1)", show_column_names = FALSE) +
      Heatmap(res$log2FoldChange[l], show_row_names = FALSE, width = unit(5, "mm"),
              name = "log2FoldChange", show_column_names = FALSE,
              col = colorRamp2(c(-2, 0, 2), c("#1984c5", "white", "#e14b31")))
    ht = draw(ht, merge_legend = TRUE)
    ht
}

# make the MA-plot with some genes highlighted
make_maplot = function(res, highlight = NULL) {
  col = rep("#00000020", nrow(res))
  cex = rep(0.5, nrow(res))
  names(col) = rownames(res)
  names(cex) = rownames(res)
  if(!is.null(highlight)) {
    col[highlight] = "red"
    cex[highlight] = 1
  }
  x = res$baseMean
  y = res$log2FoldChange
  y[y > 2] = 2
  y[y < -2] = -2
  col[col == "red" & y < 0] = "darkgreen"
  par(mar = c(4, 4, 1, 1))
  
  suppressWarnings(
    plot(x, y, col = col, 
         pch = ifelse(res$log2FoldChange > 2 | res$log2FoldChange < -2, 1, 16), 
         cex = cex, log = "x",
         xlab = "baseMean", ylab = "log2 fold change")
  )
}

# make the volcano plot with some genes highlited
make_volcano = function(res, highlight = NULL) {
  col = rep("#00000020", nrow(res))
  cex = rep(0.5, nrow(res))
  names(col) = rownames(res)
  names(cex) = rownames(res)
  if(!is.null(highlight)) {
    col[highlight] = "red"
    cex[highlight] = 1
  }
  x = res$log2FoldChange
  y = -log10(res$padj)
  col[col == "red" & x < 0] = "darkgreen"
  par(mar = c(4, 4, 1, 1))
  
  suppressWarnings(
    plot(x, y, col = col, 
         pch = 16, 
         cex = cex,
         xlab = "log2 fold change", ylab = "-log10(FDR)")
  )
  
}

#make_heatmap()
ht <- make_heatmap()
#make_maplot(res)
#make_volcano(res)
make_heatmap(filter_method=filter_by_lfc)

tiff(paste(analysis_name, 'heatmap.sig_genes.tiff', sep='.'), width=1000, height=1250)
make_heatmap()
dev.off()

library(shiny)
library(DT)
library(shinydashboard)
body = dashboardBody(
  fluidRow(
    column(width = 4,
           box(title = "Differential heatmap", width = NULL, solidHeader = TRUE, status = "primary",
               originalHeatmapOutput("ht", height = 800, containment = TRUE)
           )
    ),
    column(width = 4,
           id = "column2",
           box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
               subHeatmapOutput("ht", title = NULL, containment = TRUE)
           ),
           box(title = "Output", width = NULL, solidHeader = TRUE, status = "primary",
               HeatmapInfoOutput("ht", title = NULL)
           ),
           box(title = "Note", width = NULL, solidHeader = TRUE, status = "primary",
               htmlOutput("note")
           ),
    ),
    column(width = 4,
           box(title = "MA-plot", width = NULL, solidHeader = TRUE, status = "primary",
               plotOutput("ma_plot")
           ),
           box(title = "Volcano plot", width = NULL, solidHeader = TRUE, status = "primary",
               plotOutput("volcano_plot")
           ),
           box(title = "Result table of the selected genes", width = NULL, solidHeader = TRUE, status = "primary",
               DTOutput("res_table")
           )
    ),
    tags$style("
            .content-wrapper, .right-side {
                overflow-x: auto;
            }
            .content {
                min-width:1500px;
            }
        ")
  )
)

#library(DT)
#library(GetoptLong) # for qq() function
brush_action = function(df, input, output, session) {

  row_index = unique(unlist(df$row_index))
  selected = env$row_index[row_index]

  output[["ma_plot"]] = renderPlot({
    make_maplot(res, selected)
  })

  output[["volcano_plot"]] = renderPlot({
    make_volcano(res, selected)
  })

  output[["res_table"]] = renderDT(
    formatRound(datatable(as.matrix(res[selected, c("baseMean", "log2FoldChange", "padj")]), rownames = TRUE), columns = 1:3, digits = 3)
  )

  output[["note"]] = renderUI({
    if(!is.null(df)) {
      HTML(qq("<p>Row indices captured in <b>Output</b> only correspond to the matrix of the differential genes. To get the row indices in the original matrix, you need to perform:</p>
<pre>
l = res$padj <= @{input$fdr} &
    res$baseMean >= @{input$base_mean} &
    abs(res$log2FoldChange) >= @{input$log2fc}
l[is.na(l)] = FALSE
which(l)[row_index]
</pre>
<p>where <code>res</code> is the complete data frame from DESeq2 analysis and <code>row_index</code> is the <code>row_index</code> column captured from the code in <b>Output</b>.</p>"))
    }
  })
}

ui = dashboardPage(
  dashboardHeader(title = "DESeq2 results"),
  dashboardSidebar(
    selectInput("fdr", label = "Cutoff for FDRs:", c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05)),
    numericInput("base_mean", label = "Minimal base mean:", value = 0),
    numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 0),
    actionButton("filter", label = "Generate heatmap")
  ),
  body
)


server = function(input, output, session) {
  observeEvent(input$filter, {
    ht = make_heatmap(fdr = as.numeric(input$fdr), base_mean = input$base_mean, log2fc = input$log2fc)
    if(!is.null(ht)) {
      makeInteractiveComplexHeatmap(input, output, session, ht, "ht",
                                    brush_action = brush_action)
    } else {
      # The ID for the heatmap plot is encoded as @{heatmap_id}_heatmap, thus, it is ht_heatmap here.
      output$ht_heatmap = renderPlot({
        grid.newpage()
        grid.text("No row exists after filtering.")
      })
    }
  }, ignoreNULL = FALSE)
}

enableBookmarking(store = "url")

#shinyApp(ui, server)

# use the consensus DEG from the triplicate experiment (10-6, 12-3, 12-4).
#   sig genes that overlap between SS vs P and CI vs P
sigGenes.consensus.df <- CONSENSUS_DEG_FILE %>%
  read_csv() %>%
  dplyr::filter((SS_P_log2FoldChange > THRESHOLD_LFC & CI_P_log2FoldChange > THRESHOLD_LFC) | 
                  SS_P_log2FoldChange < -THRESHOLD_LFC & CI_P_log2FoldChange < -THRESHOLD_LFC
                ) %>%
  dplyr::filter(SS_P_padj < THRESHOLD_P & CI_P_padj < THRESHOLD_P)

sigGenes.consensus <- sigGenes.consensus.df %>%
  pull(gene_name) %>%
  unique()
sigGenes.consensus %>% head()

upGenes.consensus <- sigGenes.consensus.df %>%
  dplyr::filter(SS_P_log2FoldChange > 0 & CI_P_log2FoldChange > 0) %>%
  pull(gene_name) %>%
  unique()

downGenes.consensus <- sigGenes.consensus.df %>%
  dplyr::filter(SS_P_log2FoldChange < 0 & CI_P_log2FoldChange < 0) %>%
  pull(gene_name) %>%
  unique()

# get corresponding gene IDs
sigGenes.consensus.geneID <- data.frame(gene_name=sigGenes.consensus) %>%
  dplyr::inner_join(geneID_to_geneName, by='gene_name') %>%
  pull(gene_id)
upGenes.consensus.geneID <- data.frame(gene_name=upGenes.consensus) %>%
  dplyr::inner_join(geneID_to_geneName, by='gene_name') %>%
  pull(gene_id)
downGenes.consensus.geneID <- data.frame(gene_name=downGenes.consensus) %>%
  dplyr::inner_join(geneID_to_geneName, by='gene_name') %>%
  pull(gene_id)



# ComplexHeatMap of these select genes
filter_by_list <- function(res, ref=sigGenes.consensus.geneID) {
  rownames(res) %in% ref
}
ht_select <- make_heatmap(filter_method=filter_by_list)
tiff(paste(analysis_name, 'heatmap.consensus_deg.tiff', sep='.'), width=1000, height=1250)
print(ht_select)
dev.off()

# upregulated in Q (consensus DEG)
filter_upregulated <- function(res, ref=upGenes.consensus.geneID) {
  filter_by_list(res, ref)
}
ht_select <- make_heatmap(filter_method=filter_upregulated)
tiff(paste(analysis_name, 'heatmap.consensus_deg.up.tiff', sep='.'), width=1000, height=1250)
print(ht_select)
dev.off()

# downregulated in Q (consensus DEG)
filter_downregulated <- function(res, ref=downGenes.consensus.geneID) {
  filter_by_list(res, ref)
}
ht_select <- make_heatmap(filter_method=filter_downregulated)
tiff(paste(analysis_name, 'heatmap.consensus_deg.down.tiff', sep='.'), width=1000, height=1250)
print(ht_select)
dev.off()



# pheatmap of select genes, similar idea as above but better for showing specific gene lists/ontologies
#   (for now: using consensus DEG from previous triplicate experiment, SS vs P AND CI vs P)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
rld <- rlog(dds)
counts.lognorm <- as.data.frame(assay(rld)) %>%
  rownames_to_column(var = "gene_id") %>%
  mutate(order = 1:nrow(.)) %>%  # Use nrow(.) instead of length(.)
  dplyr::inner_join(., geneID_to_geneName, by = "gene_id") %>%
  arrange(order) %>% dplyr::select(-order)
head(counts.lognorm)

duplicates <- counts.lognorm %>%
  group_by(gene_name) %>%       # Group by 'gene_name'
  dplyr::filter(dplyr::n() > 1) %>%           # Keep only groups with more than one occurrence
  ungroup()                     # Ungroup the data for cleaner output

# View the duplicates
print(as.data.frame(duplicates) %>% head())

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
counts_norm_filter <- counts.lognorm %>% dplyr::filter(gene_name %in% sigGenes.consensus) %>%
  group_by(gene_name) %>%        # Group by gene_name
  dplyr::filter(dplyr::n() == 1) %>%           # Keep only groups with exactly one occurrence
  ungroup() %>%
  dplyr::select(-gene_id) %>%
  remove_rownames %>% column_to_rownames(var="gene_name")
as_tibble(counts_norm_filter)
nrow(counts_norm_filter)
length(sigGenes.consensus)

counts_norm_filter_zscore <- t(apply(counts_norm_filter, 1, cal_z_score))
sample_col <- data.frame(sample=samples[[compare_by$col]])
pheatmap.out <- pheatmap(counts_norm_filter_zscore, cluster_cols=F, annotation_col=sample_col) # TODO: clean this up
pheatmap.out
tiff(paste(analysis_name, 'pheatmap.consensus_deg.tiff', sep='.'), width=1000, height=1250)
print(pheatmap.out)
dev.off()

# chromatin remodelers
chromatin_remodelers.file <- CHROMATIN_REMODELERS_FILE
chromatin_remodelers <- read.table(chromatin_remodelers.file) %>% pull(V1)
chromatin_remodelers %>% head()

counts_norm_filter <- counts.lognorm %>% 
  dplyr::filter(gene_name %in% sigGenes.consensus & gene_name %in% chromatin_remodelers) %>%
  group_by(gene_name) %>%        # Group by gene_name
  dplyr::filter(dplyr::n() == 1) %>%           # Keep only groups with exactly one occurrence
  ungroup() %>%
  dplyr::select(-gene_id) %>%
  remove_rownames %>% column_to_rownames(var="gene_name")
as_tibble(counts_norm_filter)
nrow(counts_norm_filter)
length(sigGenes.consensus)

counts_norm_filter_zscore <- t(apply(counts_norm_filter, 1, cal_z_score))
sample_col <- data.frame(sample=samples[[compare_by$col]])
# pheatmap.out <- pheatmap(counts_norm_filter_zscore, cluster_cols=F, annotation_col=sample_col) # TODO: clean this up
# pheatmap.out
# 
# tiff(paste(analysis_name, 'pheatmap.consensus_deg.chromatin_remodelers.tiff', sep='.'), width=1000, height=1250)
# print(pheatmap.out)
# dev.off()

# genes that are significantly up/downregulated
sig_features <- na.omit(res)
downregulated_features <- rownames(sig_features[sig_features$log2FoldChange < -THRESHOLD_LFC,]) %>%
  sub("\\.[0-9]+$", "", .) %>% unique()
upregulated_features <- rownames(sig_features[sig_features$log2FoldChange > THRESHOLD_LFC,]) %>%
  sub("\\.[0-9]+$", "", .) %>% unique()
sig_features <- rownames(sig_features[abs(sig_features$log2FoldChange) > THRESHOLD_LFC,]) %>%
  sub("\\.[0-9]+$", "", .) %>% unique()
universe_features <- rownames(res) %>%
  sub("\\.[0-9]+$", "", .) %>% unique()

# compare to consensus DEG list
sigGenes.consensus.df
res.compare_to_consensus <- res[rownames(res) %in% unique(sigGenes.consensus.df$gene_id.x),]
length(res.compare_to_consensus$log2FoldChange)
length(sigGenes.consensus)



# GENE ONTOLOGY----------
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# GO analysis
# TODO: try out PANTHER as well and see which is better suited
GO_upregulated <- enrichGO(gene = upregulated_features, universe = universe_features, OrgDb="org.Hs.eg.db", keyType="ENSEMBL", ont="ALL")
print(as_tibble(GO_upregulated), n=50)
fit <- barplot(GO_upregulated, showCategory=20, font=10)
plot(fit)
ggsave(paste(analysis_name, 'GO.up.tiff', sep='.'), plot = fit, width = 10, height = 10, dpi = 300)

GO_downregulated <- enrichGO(gene = downregulated_features, universe = universe_features, OrgDb="org.Hs.eg.db", keyType="ENSEMBL", ont="ALL")
print(as_tibble(GO_downregulated), n=50)
fit <- barplot(GO_downregulated, showCategory=20, font=10)
plot(fit)
ggsave(paste(analysis_name, 'GO.down.tiff', sep='.'), plot = fit, width = 10, height = 10, dpi = 300)

GO_sig <- enrichGO(gene = sig_features, universe = universe_features, OrgDb="org.Hs.eg.db", keyType="ENSEMBL", ont="ALL")
print(as_tibble(GO_sig), n=50)
fit <- barplot(GO_sig, showCategory=20, font=10)
plot(fit)
ggsave(paste(analysis_name, 'GO.sig.tiff', sep='.'), plot = fit, width = 10, height = 10, dpi = 300)



# optional: run heatmaps.R to identify gene expression patterns of genes differentially expressed between proliferation and quiescence