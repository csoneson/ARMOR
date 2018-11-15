args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## This script performs differential transcript usage analysis with DRIMSeq,
## based on abundance estimates from Salmon.
## Please modify the design, contrast, filtering criterion and the desired level
## (gene or feature) below as appropriate

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(tximeta))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DRIMSeq))
suppressPackageStartupMessages(library(ggplot2))

## Open pdf file to contain any figure generated below
pdf(gsub("rds$", "pdf", outrds))

## Load the SummarizedExperiment object obtained from tximeta
se <- readRDS(se)

## Quantification on the gene level
sg <- se$sg

## Quantification on the transcript level
st <- se$st

## Create dmDSdata object
counts <- data.frame(feature_id = rownames(st),
                     gene_id = unlist(rowData(st)$gene_id),
                     assays(st)[["counts"]],
                     row.names = NULL)

metadata <- metadata %>% select(sample_id = names, group  = celline)
  
d <- dmDSdata(counts = counts, samples = metadata)
plotData(d)

## Filter ************** MODIFY ************** 
d <- dmFilter(d, min_samps_gene_expr = 3, min_samps_feature_expr = 3,
              min_gene_expr = 10, min_feature_expr = 5)

## Define design. ************** MODIFY ************** 
# (des <- model.matrix(~ XXXX, data = samples(d)))
(des <- model.matrix(~ 0 + group, data = samples(d)))
## Calculate precision
set.seed(123)
d <- dmPrecision(d, design = des)
plotPrecision(d)

## Fit model
d <- dmFit(d, design = des, verbose = 1)

## Define contrasts. ************** MODIFY ************** 
#(contrasts <- as.data.frame(makeContrasts(XXXX, levels = des)))
(contrasts <- as.data.frame(makeContrasts("groupN61311-groupN052611", levels = des)))
## Perform tests ************** MODIFY **************
level <- "gene" ## set to "feature" if transcript-level results are desired
signif3 <- function(x) signif(x, digits = 3)
DRIMSeq_fits <- lapply(contrasts, function(cm) {
  dr <- dmTest(d, contrast = cm, verbose = 1)
  print(plotPValues(dr, level = level))
  dr
})

DRIMSeq_res <- lapply(DRIMSeq_fits, function(dr) {
  results(dr, level = level) %>% dplyr::mutate_if(is.numeric, signif3)
})

## Write results to text files
if (class(DRIMSeq_res) == "data.frame") {
  write.table(DRIMSeq_res %>% dplyr::arrange(pvalue), 
              file = gsub("rds$", "txt", outrds), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
} else {
  for (nm in names(DRIMSeq_res)) {
    write.table(DRIMSeq_res[[nm]] %>% dplyr::arrange(pvalue), 
                file = gsub("\\.rds$", paste0("_", nm, ".txt"), outrds), 
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}

dev.off()

saveRDS(list(results = DRIMSeq_res, fits = DRIMSeq_fits), file = outrds)

sessionInfo()
date()
