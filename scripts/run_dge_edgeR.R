args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## This script performs differential expression analysis with edgeR, based on 
## abundance estimates from Salmon. It supports testing one or more contrasts. 
## To be compatible with downstream analysis scripts, the output has to be
## either:
## - a data frame with results (e.g., returned by edgeR's topTags(...)$table).
## There must be one column named "gene" that gives the gene ID.
## - a named list of such data frames, one for each contrast (recommended).
## To run the script, modify at least the definition of the design matrix and
## the contrasts of interest.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(tximeta))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(ggplot2))

print(se)

## Open pdf file to contain any figure generated below
pdf(gsub("rds$", "pdf", outrds))

## Load the SummarizedExperiment object obtained from tximeta
se <- readRDS(se)

## Use the quantification on the gene level
sg <- se$sg
metadata <- colData(sg)

## Create DGEList and include average transcript length offsets
cts <- assays(sg)[["counts"]]
normMat <- assays(sg)[["length"]]
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
dge0 <- DGEList(cts)
dge0$offset <- t(t(log(normMat)) + o)
dge0 <- calcNormFactors(dge0)

## Add gene annotation to the DGEList object
dge0$genes <- as.data.frame(rowRanges(sg))

## Define design. ************** MODIFY ************** 
stopifnot(all(colnames(dge0) == metadata$names))
#(des <- model.matrix(~ XXXX, data = metadata))
#e.g.
(des <- model.matrix(~ 0 + celline, data = metadata))

## Filter out genes with average CPM below 1
print(dim(dge0))
cpms <- cpm(dge0)
dge <- dge0[apply(cpms, 1, mean) > 1, ]
dge <- calcNormFactors(dge)
print(dim(dge))

## Estimate dispersion and fit model
dge <- estimateDisp(dge, design = des)
qlfit <- glmQLFit(dge, design = des)

## Plot dispersions
plotBCV(dge)

## Define contrasts. ************** MODIFY ************** 
#(contrasts <- as.data.frame(makeContrasts(XXXX, levels = des)))
#e.g.
(contrasts <- as.data.frame(makeContrasts(contrasts= "cellineN61311-cellineN052611", levels = des)))

## Perform tests
signif3 <- function(x) signif(x, digits = 3)
edgeR_res <- lapply(contrasts, function(cm) {
  qlf <- glmQLFTest(qlfit, contrast = cm)
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  tt %>% dplyr::mutate_if(is.numeric, signif3)
})

## Write results to text files and make MA plots
if (class(edgeR_res) == "data.frame") {
  write.table(edgeR_res %>% dplyr::arrange(PValue), 
              file = gsub("rds$", "txt", outrds), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  print(ggplot(edgeR_res, aes(x = logCPM, y = logFC, color = FDR <= 0.05)) + 
          geom_point() + theme_bw() + 
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")))
} else {
  for (nm in names(edgeR_res)) {
    write.table(edgeR_res[[nm]] %>% dplyr::arrange(PValue), 
                file = gsub("\\.rds$", paste0("_", nm, ".txt"), outrds), 
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    print(ggplot(edgeR_res[[nm]], aes(x = logCPM, y = logFC, color = FDR <= 0.05)) + 
            geom_point() + theme_bw() + 
            scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + 
            ggtitle(nm))
  }
}

dev.off()

saveRDS(list(results = edgeR_res, data = dge0), file = outrds)

sessionInfo()
date()
