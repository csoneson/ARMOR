args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tximport)
  library(tximeta)
  library(SummarizedExperiment)
})

print(salmondir)
print(json)
print(metafile)
print(outrds)

## Load json linkedTxome
loadLinkedTxome(json)

## Read metadata
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")

## List Salmon directories
salmonfiles <- paste0(salmondir, "/", metadata$names, "/quant.sf")
names(salmonfiles) <- metadata$names

## Add file column to metadata and import annotated abundances
## In transcript level
coldata <- cbind(metadata, files = salmonfiles, stringsAsFactors = FALSE)
st <- tximeta(coldata)

## Summarize to gene level
sg <- summarizeToGene(st)

## Add gene information, e.g. gene_name, entrezid, ... to se
rowData(st) <- rowData(st) %>%
    data.frame() %>%
    left_join(data.frame(rowData(sg))) %>%
    DataFrame()

saveRDS(list(st = st, sg = sg), file = outrds)

sessionInfo()
date()


