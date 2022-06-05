args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(dplyr)
    library(tximport)
    library(tximeta)
    library(SingleCellExperiment)
})

print(salmondir)
print(json)
print(metafile)
print(outrds)
print(annotation)
print(organism)

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
st <- tximeta::tximeta(coldata)

## Summarize to gene level
sg <- summarizeToGene(st)

## If the 'entrezid' column exists and is a list, convert to a vector
if ("entrezid" %in% colnames(rowData(sg)) && 
    is(rowData(sg)$entrezid, "list")) {
    if (any(vapply(rowData(sg)$entrezid, length, 1) > 1)) {
        warning("Some elements of rowData(sg)$entrezid consisted of ",
                "more than one object. Only the first one is retained.")
    }
    rowData(sg)$entrezid <- vapply(
        rowData(sg)$entrezid, 
        function(w) w[[1]], 
        as(NA, class(rowData(sg)$entrezid[[1]]))
    )
}

## Add gene_names for Gencode reference
if(annotation == "Gencode") {
    if(organism == "Homo_sapiens") {
        library(org.Hs.eg.db)
    } else {
        library(org.Mm.eg.db)
    }
    sg <- tximeta::addIds(sg, "SYMBOL", gene = TRUE)
    rowData(sg)$gene_name <- rowData(sg)$SYMBOL
}

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
    dplyr::left_join(data.frame(rowData(sg))) %>%
    DataFrame()

## Change the row names in sg to have geneID__geneSymbol
rownames(sg) <- paste(rowData(sg)$gene_id, rowData(sg)$gene_name, sep = "__")

# Coerce the object from SummarizedExperiment to SingleCellExperiment
st <- as(st, "SingleCellExperiment")
sg <- as(sg, "SingleCellExperiment")

saveRDS(list(st = st, sg = sg), file = outrds)

sessionInfo()
date()


