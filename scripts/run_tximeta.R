args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
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

## If rowData(st)$gene_id is a CharacterList, convert it to character to allow 
## the joining below
if (is(rowData(st)$gene_id, "CharacterList")) {
    if (any(vapply(rowData(st)$gene_id, length, 1) > 1)) {
        warning("Some elements of rowData(st)$gene_id consisted of more than one",
                "object. Only the first one is retained.")
    }
    rowData(st)$gene_id <- vapply(rowData(st)$gene_id, function(w) w[[1]], "")
}

## If rowData(st)$tx_id is of class integer, replace it with the tx_name 
## column
if (is(rowData(st)$tx_id, "integer")) {
    rowData(st)$tx_id <- rowData(st)$tx_name
}

## Add gene information, e.g. gene_name, entrezid, ... (if provided) to
## transcript-level SE
rowData(st) <- rowData(st) %>%
    data.frame() %>%
    left_join(data.frame(rowData(sg))) %>%
    DataFrame()

saveRDS(list(st = st, sg = sg), file = outrds)

sessionInfo()
date()


