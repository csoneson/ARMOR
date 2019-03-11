args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(tximeta)
})

print(transcriptfasta)
print(salmonidx)
print(gtf)

print(annotation)
ss <- strsplit(organism, "_")[[1]]
organism <- paste(paste(ss[1], ss[2]))
print(organism)
print(release)
print(build)
print(output)

makeLinkedTxome(indexDir = dirname(salmonidx),
                source = annotation,
                organism = organism,
                release = release,
                genome = build,
                fasta = transcriptfasta,
                gtf = gtf,
                write = TRUE, 
                jsonFile = output)

sessionInfo()
date()
