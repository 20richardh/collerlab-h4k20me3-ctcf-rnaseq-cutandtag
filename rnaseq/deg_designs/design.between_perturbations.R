### perturbations such as A196 vs DMSO (control), CTCF_OE vs RFP_OE, siCTCF vs siControl
### CTCF_OE supposedly makes cells more proliferating. We want to test this
###   vs the RFP_OE control (OE protocol does seem to cause some outlier outcomes
###   compared to the other conditions, so we're just doing OE vs OE)
### ADDITIONALLY, because CTCF_OE,CI runs didn't really overexpress, we are
###   going to just look at SS/P juxtaposition (OE worked well)

SKIP_PCA <- F
# REF <- 'RFP_OE'
# PERTURB <- 'CTCF_OE'
# REF <- 'SiControl'
# PERTURB <- 'SiCTCF'
REF <- 'DMSO'
PERTURB <- 'A196'


compare_by <- data.frame(col=c('gene_expression'), ref=c(REF), perturb=c(PERTURB))

# compare_by <- data.frame(col=c('bioReplicate_id'), ref=c('12_3'), perturb=c('12_4'))
# sort_col <- 'bioReplicate_id'

comparison_name <- paste(compare_by$perturb, compare_by$ref, sep='_vs_')
comparison_name

### by gene expression feature
# col_of_interest <- 'quiescence'
#features_of_interest <- c('P', 'Q') ### see comments at top for OE experiments
col_of_interest <- 'cell_type'
features_of_interest <- c('P', 'SS')
# features_of_interest <- c('SS')
sort_col <- 'gene_expression_idx'
to_exclude <- list()