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

print(salmondir)
print(json)
print(metafile)
print(outrds)

## Load json linkedTxome
loadLinkedTxome(json)

## Open pdf file to contain any figure generated below
pdf(gsub("rds$", "pdf", outrds))

## Read metadata
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")

## List Salmon directories
salmonfiles <- paste0(salmondir,"/",metadata$names, "/quant.sf")
names(salmonfiles) <- metadata$names
(salmonfiles <- salmonfiles[file.exists(salmonfiles)])

## Add file column to metadata and import annotated abundances
coldata <- cbind(metadata, files = salmonfiles, stringsAsFactors=FALSE)
se <- tximeta(coldata)

## Create dmDSdata object
counts <- data.frame(feature_id = rownames(se),
                     gene_id = unlist(rowData(se)$gene_id),
                     assays(se)[["counts"]],
                     row.names = NULL)

metadata <- metadata %>% select(sample_id = names, group  = celline)
  
d <- dmDSdata(counts = counts, samples = metadata)
plotData(d)

## Filter ************** MODIFY ************** 
d <- dmFilter(d, min_samps_gene_expr = 3, min_samps_feature_expr = 3,
              min_gene_expr = 10, min_feature_expr = 5)

## Define design. ************** MODIFY ************** 
(des <- model.matrix(~ XXXX, data = samples(d)))

## Calculate precision
set.seed(123)
d <- dmPrecision(d, design = des)
plotPrecision(d)

## Fit model
d <- dmFit(d, design = des, verbose = 1)

## Define contrasts. ************** MODIFY ************** 
(contrasts <- as.data.frame(makeContrasts(XXXX, levels = des)))

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
