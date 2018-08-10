args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages(library(tximeta)) 

print(transcriptfasta)
print(salmonidx)
print(gtf)

print(annotation)
print(organism)
ss <- strsplit(organism, "_")[[1]]
organism <- paste(paste(ss[1], ss[2]))
print(release)
print(build)

makeLinkedTxome(indexDir=salmonidx,
                source=annotation,
                organism=organism,
                release=release,
                genome=build,
                fasta=transcriptfasta,
                gtf=gtf,
                write=FALSE)

sessionInfo()
date()
