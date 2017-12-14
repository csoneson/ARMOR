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
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DRIMSeq))
suppressPackageStartupMessages(library(ggplot2))

print(tx2gene)
print(salmondir)
print(metafile)
print(outrds)

## Open pdf file to contain any figure generated below
pdf(gsub("rds$", "pdf", outrds))

## List Salmon directories
salmondirs <- list.files(salmondir, full.names = TRUE)
salmonfiles <- paste0(salmondirs, "/quant.sf")
names(salmonfiles) <- basename(salmondirs)
(salmonfiles <- salmonfiles[file.exists(salmonfiles)])

## Read transcript-to-gene mapping
tx2gene <- readRDS(tx2gene)

## Read Salmon abundances
txi <- tximport(files = salmonfiles, type = "salmon", txOut = TRUE, 
                dropInfReps = TRUE)
counts <- as.data.frame(txi$counts) %>% 
  tibble::rownames_to_column("feature_id") %>%
  dplyr::left_join(tx2gene %>% dplyr::select(gene, tx), by = c("feature_id" = "tx")) %>%
  dplyr::rename(gene_id = gene) %>%
  dplyr::select(gene_id, feature_id, everything())

## Read metadata and reorder in the same order as the count matrix
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")
stopifnot(all(metadata$ID %in% colnames(counts)))
stopifnot(all(setdiff(colnames(counts), c("gene_id", "feature_id")) %in% metadata$ID))
metadata <- metadata[match(setdiff(colnames(counts), c("gene_id", "feature_id")), 
                           metadata$ID), ] %>%
  dplyr::rename(sample_id = ID)

## Create dmDSdata object
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
