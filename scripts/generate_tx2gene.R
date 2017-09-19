args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## Parse a transcript fasta file from Ensembl and extract information about the
## transcripts.

suppressPackageStartupMessages(library(Biostrings))

print(transcriptfasta)
print(outrds)

transcripts <- readDNAStringSet(transcriptfasta)

tx2gene <- data.frame(t(sapply(as.character(names(transcripts)), function(nm) {
  a <- strsplit(nm, " ")[[1]]
  tx <- a[1]
  gene <- gsub("^gene:", "", a[grep("^gene:", a)])
  symbol <- gsub("^gene_symbol:", "", a[grep("^gene_symbol:", a)])
  gene_biotype <- gsub("^gene_biotype:", "", a[grep("^gene_biotype:", a)])
  tx_biotype <- gsub("^transcript_biotype:", "", a[grep("^transcript_biotype:", a)])
  position <- gsub("chromosome:", "", a[grep("^chromosome:", a)])
  c(tx = ifelse(length(tx) != 0, tx, NA), 
    gene = ifelse(length(gene) != 0, gene, NA),
    symbol = ifelse(length(symbol) != 0, symbol, NA),
    gene_biotype = ifelse(length(gene_biotype) != 0, gene_biotype, NA),
    tx_biotype = ifelse(length(tx_biotype) != 0, tx_biotype, NA),
    chromosome = ifelse(length(position) != 0, strsplit(position, ":")[[1]][2], NA),
    start = ifelse(length(position) != 0, strsplit(position, ":")[[1]][3], NA),
    end = ifelse(length(position) != 0, strsplit(position, ":")[[1]][4], NA),
    strand = ifelse(length(position) != 0, strsplit(position, ":")[[1]][5], NA))
})), stringsAsFactors = FALSE)
rownames(tx2gene) <- NULL

saveRDS(tx2gene, file = outrds)

sessionInfo()
date()
