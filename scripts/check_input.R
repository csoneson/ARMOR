args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(edgeR)
})

## ----------------The input arguments------------------------------------------
print(metafile)
print(outFile)
print(gtf)
print(genonme)
print(fastqdir)
print(fqsuffix)

msg0 <- try({
    pt <- paste0("\\.", fqsuffix, ".1gz")
    lf <- list.files(pattern = pt, path = fastqdir)
    if (length(lf) == 0) {
        stop("FASTQ files with ", fqsuffix, " .gz are not available. \n")
    }
}, silent = TRUE)

msg1 <- try({
    fe <- file.exists(genonme)
    if (!fe) {
        stop("The genonme file doesn't exist. \n")
    }
}, silent = TRUE)

msg2 <- try({
    fe <- file.exists(gtf)
    if (!fe) {
        stop("The gtf file doesn't exist. \n")
    }
}, silent = TRUE)

msg3 <- try({
    if (exists("design")) {
        print(design)
    } else {
        design <- NULL
    }
}, silent = TRUE)


msg4 <- try({
    if (exists("contrast")) {
        contrast <- strsplit(gsub(" ","",contrast), ",")[[1]]
        print(contrast)
    } else {
        contrast <- NULL
    }
}, silent = TRUE)

## ---------------------------Test run -------------------------------
## Read metadata
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")
rownames(metadata) <- metadata$names
metadata

## Define design matrix
msg5 <- try({des <- model.matrix(as.formula(design), data = metadata)},
            silent = TRUE)

# Define contrasts
msg6 <- try({contrasts <- as.data.frame(makeContrasts(contrasts = contrast, 
                                                      levels = des))},
            silent = TRUE)

msgL <- list(msg0, msg1, msg2, msg3, msg4, msg5, msg6)
isError <- lapply(msgL, FUN = function(x) {class(x) == "try-error"})
isError <- unlist(isError)

msg <- msgL[isError]

if (length(msg) > 0) {
    capture.output(writeLines(msg[[1]]), file = outFile)
} else {
    mylist <- list("Design matrix" = des, "Contrasts matrix" = contrasts)
    capture.output(mylist, file = outFile)
}



