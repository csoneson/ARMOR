args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages(library(tximeta)) #how to install this?
suppressPackageStartupMessages(library(SummarizedExperiment))

print(transcriptfasta)
print(outrds)
print(annotation)
print(organism)
print(release)
print(build)
print(salmondir)
print(salmonidx)
print(gtf)

salmondirs <- list.files(salmondir, full.names = TRUE)
salmonfiles <- paste0(salmondirs, "/quant.sf")

coldata <- data.frame(files = salmonfiles, 
                      names = basename(salmondirs), 
                      stringsAsFactors=FALSE)

se <- tximeta(coldata)

if(ncol(rowData(se)) == 0){
  makeLinkedTxome(indexDir=salmonidx,
                  source=annotation,
                  organism=organism,
                  release=release,
                  genome=build,
                  fasta=transcriptfasta,
                  gtf=gtf,
                  write=FALSE)
  se <- tximeta(coldata)
}

tx2gene <- data.frame(tx = rowData(se)$tx_id,
                      gene = rowData(se)$gene_id,
                      #symbol = mcols(rowRanges(se))$symbol,
                      gene_biotype = rowData(se)$tx_biotype,
                      tx_biotype = rowData(se)$tx_biotype, #cheat
                      chromosome = as.character(seqnames(se)),
                      start = start(se),
                      end = end(se),
                      strand = as.character(strand(se)),
                      stringsAsFactors = FALSE)

saveRDS(tx2gene, file = outrds)

sessionInfo()
date()
