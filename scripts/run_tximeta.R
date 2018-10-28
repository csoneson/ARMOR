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

print(salmondir)
print(json)
print(metafile)
print(outrds)

## Load json linkedTxome
loadLinkedTxome(json)

## Read metadata
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")

## List Salmon directories
salmonfiles <- paste0(salmondir,"/",metadata$names, "/quant.sf")
names(salmonfiles) <- metadata$names
(salmonfiles <- salmonfiles[file.exists(salmonfiles)])

## Add file column to metadata and import annotated abundances
coldata <- cbind(metadata, files = salmonfiles, stringsAsFactors=FALSE)
se <- tximeta(coldata)

## Summarize to gene level
sg <- summarizeToGene(se)
