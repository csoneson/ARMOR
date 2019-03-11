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

msg0 <- try({
    if (exists("design")) {
        print(design)
    } else {
        design <- NULL
    }
}, silent = TRUE)


msg1 <- try({
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
msg2 <- try({des <- model.matrix(as.formula(design), data = metadata)},
            silent = TRUE)

# Define contrasts
msg3 <- try({contrasts <- as.data.frame(makeContrasts(contrasts = contrast, 
                                                      levels = des))},
            silent = TRUE)

msgL <- list(msg0, msg1, msg2, msg3)
isError <- lapply(msgL, FUN = function(x) {class(x) == "try-error"})
isError <- unlist(isError)

msg <- msgL[isError]

if (length(msg) > 0) {
    capture.output(writeLines(msg[[1]]), file = outFile)
} else {
    mylist <- list("Design matrix" = des, "Contrasts matrix" = contrasts)
    capture.output(mylist, file = outFile)
}



