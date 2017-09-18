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
## - a named list of such data frames, one for each contrast.
## To run the script, please modify at least the definition of the design matrix
## and the contrasts of interest.

suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(edgeR))

print(tx2gene)
print(salmondir)
print(metafile)
print(outrds)

## List Salmon directories
salmondirs <- list.files(salmondir, full.names = TRUE)
salmonfiles <- paste0(salmondirs, "/quant.sf")
names(salmonfiles) <- basename(salmondirs)
(salmonfiles <- salmonfiles[file.exists(salmonfiles)])

## Read transcript-to-gene mapping
tx2gene <- readRDS(tx2gene)

## Read Salmon abundances
txi <- tximport(files = salmonfiles, type = "salmon", txOut = FALSE, 
                tx2gene = tx2gene[, c("tx", "gene")])

## Read metadata and reorder in the same order as count matrix
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")
metadata <- metadata[match(colnames(txi$counts), metadata$ID), ]

## Create DGEList and include average transcript length offsets
cts <- txi$counts
normMat <- txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
dge0 <- DGEList(cts)
dge0$offset <- t(t(log(normMat)) + o)
dge0 <- calcNormFactors(dge0)

## Define design. MODIFY
stopifnot(all(colnames(dge0) == metadata$ID))
des <- model.matrix(~ XXXX, data = metadata)

## Filter out genes with average CPM below 1
print(dim(dge0))
cpms <- cpm(dge0)
dge <- dge0[apply(cpms, 1, mean) > 1, ]
dge <- calcNormFactors(dge)
print(dim(dge))

## Add gene annotation
annot <- tx2gene[match(rownames(dge), tx2gene$gene), ]
rownames(annot) <- annot$gene
dge$genes <- annot

## Estimate dispersion and fit model
dge <- estimateDisp(dge, design = des)
qlfit <- glmQLFit(dge, design = des)

## Define contrasts. MODIFY
(contrasts <- as.data.frame(makeContrasts(XXXX, levels = des)))

## Perform tests
edgeR_res <- lapply(contrasts, function(cm) {
  qlf <- glmQLFTest(qlfit, contrast = cm)
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  tt <- signif(tt, digits = 3)
  tt$gene <- rownames(tt)
  tt
})

saveRDS(list(results = edgeR_res, data = dge0), file = outrds)


sessionInfo()
date()
