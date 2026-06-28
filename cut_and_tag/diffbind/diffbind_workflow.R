# DiffBind differential peak analysis — Figures 2C and 2D.
# Set PROJECT_DIR in config.R before running.
#-------------------------------------------------------------------------------
source(file.path(dirname(sys.frame(1)$ofile), "config.R"))
library(GenomicRanges)

setup <- function() {
  setwd(PROJECT_DIR)
  setwd("diffbind")
  setwd("final")
}


# P,Q CTCF----------------------------------------------------------------------
SAMPLESHEET <- SAMPLESHEET_QP
setup()
COMPARISON_NAME <- "P.Q"
FACTORS <- c("CTCF")
SELECT_BY <- 'Condition'
SELECT_VALS <- c('P', 'Q')
subdir <- file.path(COMPARISON_NAME, paste(FACTORS, collapse='.'))
dir.create(COMPARISON_NAME, showWarnings = FALSE)
# contrast design for DiffBind
get_contrast <- function() {
  DBA_CONDITION
}
get_consensus <- function() {
  c(DBA_TISSUE) # c(DBA_TISSUE, DBA_CONDITION)
}
make_contrast <- function(dbObj) {
  dba.contrast(dbObj, categories=DBA_CONDITION,
               minMembers = 2)
}

dir.create(subdir, showWarnings = FALSE)
setwd(subdir)

#-------------------------------------------------------------------------------
# Differential peak analysis using DiffBind + DESeq2 workflow
#-------------------------------------------------------------------------------
library(tidyverse)
library(DiffBind)


# READ SAMPLE DESIGN
samples <- read.csv(SAMPLESHEET)
samples <- samples[samples$Factor %in% FACTORS, ]
if (!is.null(SELECT_BY)) {
  samples <- samples[samples[[SELECT_BY]] %in% SELECT_VALS,]
}
as_tibble(samples)
dbObj <- dba(sampleSheet=samples)
dbObj


# ANALYSIS----------------------------------------------------------------------

dba_consensus <- dba.peakset(dbObj, consensus=get_consensus(), minOverlap=2) # minOverlap=0.66, reevaluate when we add different peak callers
dba_consensus <- dba(dba_consensus, mask=dba_consensus$masks$Consensus, minOverlap=1)
dba_consensus
peaks_consensus <- dba.peakset(dba_consensus, bRetrieve=TRUE)
peaks_consensus

dbObj <- dba(sampleSheet=samples) %>%
  dba.blacklist() %>%
  dba.count(peaks=peaks_consensus) %>%
  dba.normalize()
dbObj <- make_contrast(dbObj) %>%
  dba.analyze(method=DBA_ALL_METHODS)
# summary of each tool, by number of differentially enriched peaks
dba_summary <- dba.show(dbObj, bContrasts=T)	
dba_summary
# which sites are differentially bound, as a report
dba.DB <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
dba.DB

normCounts <- dba.peakset(dbObj, bRetrieve=TRUE)
normCounts

# WRITE ENRICHMENT ANALYSIS
# consensus peaks
write.table(peaks_consensus %>% as.data.frame() %>% dplyr::select(seqnames, start, end), file=paste(COMPARISON_NAME, "consensus.bed", sep='.'), sep='\t', quote=F, row.names=F, col.names=F)
# deseq2 data
out <- as.data.frame(dba.DB)
write.table(out, file=paste(COMPARISON_NAME,"deseq2.txt", sep='.'), sep="\t", quote=F, row.names=F)
# enriched in set a (a vs b comparison)
enrich_a <- out %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  dplyr::select(seqnames, start, end)
write.table(enrich_a, file=paste(COMPARISON_NAME, dba_summary$Group, "bed", sep='.'), sep="\t", quote=F, row.names=F, col.names=F)
# enriched in set b (a vs b comparison)
enrich_b <- out %>% 
  filter(FDR < 0.05 & Fold < 0) %>% 
  dplyr::select(seqnames, start, end)
write.table(enrich_b, file=paste(COMPARISON_NAME, dba_summary$Group2, "bed", sep='.'), sep="\t", quote=F, row.names=F, col.names=F)


# Venn diagram
if (length(enrich_a[[1]]) > 0 && length(enrich_b[[1]]) > 0) {
  venn <- dba.plotVenn(dbObj, contrast=1, bDB=TRUE,
               bGain=TRUE, bLoss=TRUE, bAll=FALSE)
  png("venn_diagram.png", height=800, width=800)
  dba.plotVenn(dbObj, contrast=1, bDB=TRUE,
               bGain=TRUE, bLoss=TRUE, bAll=FALSE)
  dev.off()
}

# MA AND VOLCANO PLOTS
dba.plotMA(dbObj, method=DBA_DESEQ2)
dba.plotMA(dbObj, bXY=TRUE)
png("maplot.deseq2.png", height=800, width=800)
dba.plotMA(dbObj, method=DBA_DESEQ2)
dev.off()
png("maplot.xy.png", height=800, width=800)
dba.plotMA(dbObj, bXY=TRUE)
dev.off()

dba.plotVolcano(dbObj)
png("volcano.png", height=800, width=800)
dba.plotVolcano(dbObj)
dev.off()

# BOXPLOTS
pvals <- tryCatch(dba.plotBox(dbObj), error=function(e) NULL)
png("boxplots.png", height=800, width=800)
tryCatch(dba.plotBox(dbObj), error=function(e) NULL)
dev.off()
if (length(enrich_a[[1]]) > 0 && length(enrich_b[[1]]) > 0) { # nonzero enrichments in both directions (gain, loss)
  pvals
  write.table(as.data.frame(pvals), file=paste(COMPARISON_NAME, "boxplots.pvals.txt", sep='.'), sep="\t")
}
# CORR
corvals <- dba.plotHeatmap(dbObj)
png("corr.png", height=800, width=800)
dba.plotHeatmap(dbObj)
dev.off()
corvals
write.table(as.data.frame(corvals), file=paste(COMPARISON_NAME, "corvals.txt", sep='.'), sep="\t")


# fold change (raw)
sum(dba.DB$Fold<0)
sum(dba.DB$Fold>0)

# HEATMAP
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
readscores <- dba.plotHeatmap(dbObj, contrast=1, correlations=FALSE,
                              scale="row", colScheme = hmap, 
                              sites=dba.DB.down)
tiff("heatmap.tiff", height=800, width=800)
dba.plotHeatmap(dbObj, contrast=1, correlations=FALSE,
                scale="row", colScheme = hmap)
dev.off()
readscores
write.table(as.data.frame(readscores), file=paste(COMPARISON_NAME,"readscores.txt", sep='.'), sep="\t", quote=F, row.names=F)


# PROFILE PLOT (HEATMAP) OF GAINED AND LOST PEAKS
# # merged by replicate
profiles <- dba.plotProfile(dbObj, sites=dba.DB.down)
plt.out <- dba.plotProfile(profiles)
tiff("profiles.merged_by_replicate.downregulated.tiff", height=800, width=500)
print(plt.out)
dev.off()

profiles <- dba.plotProfile(dbObj, sites=peaks_consensus)
dba.plotProfile(profiles)
plt.out
tiff("profiles.merged_by_replicate.consensus.tiff", height=800, width=500)
print(plt.out)
dev.off()

# # no merging
profiles <- dba.plotProfile(dbObj, merge=NULL, sites=dba.DB.down)
dba.plotProfile(profiles)
tiff("profiles.downregulated_only.sorted.tiff", height=800, width=800)
dba.plotProfile(profiles)
dev.off()
