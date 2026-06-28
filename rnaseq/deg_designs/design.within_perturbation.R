#run_idx <- 1
compare_by <- data.frame(col=c('cell_type'), ref=c('P'), perturb=c('SS'))
#compare_by <- data.frame(col=c('quiescence'), ref=c('P'), perturb=c('Q'))

# compare_by <- data.frame(col=c('bioReplicate_id'), ref=c('12_3'), perturb=c('12_4'))
# sort_col <- 'bioReplicate_id'

comparison_name <- paste(compare_by$perturb, compare_by$ref, sep='_vs_')
comparison_name

### by gene expression feature
col_of_interest <- 'gene_expression'
features_of_interest <- "SiControl" #c(names(GENE_EXPR_COLS)[run_idx])
sort_col <- 'quiescence_idx'
to_exclude <- list()# list(list(col='gene_expression', features=c('RFP_OE', 'CTCF_OE')))