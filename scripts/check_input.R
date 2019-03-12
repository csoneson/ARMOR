
args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(edgeR)
})

## ----------------The input arguments------------------------------------------
if (exists("metafile")) {
    print(metafile) 
    } else {
        metafile <- NULL
    }

if (exists("outFile")) {
    print(outFile) 
    } else {
        outFile <- NULL
    }

if (exists("gtf")) {
    print(gtf) 
} else {
    gtf <- NULL
}

if (exists("gtf")) {
    print(gtf) 
} else {
    gtf <- NULL
}

if (exists("genome")) {
    print(genome) 
} else {
    genome <- NULL
}

if (exists("fastqdir")) {
    print(fastqdir) 
} else {
    fastqdir <- NULL
}

if (exists("fqsuffix")) {
    print(fqsuffix) 
} else {
    fqsuffix <- NULL
}

if (exists("fqext1")) {
    print(fqext1) 
} else {
    fqext1 <- NULL
}

if (exists("fqext2")) {
    print(fqext2) 
} else {
    fqext2 <- NULL
}
if (exists("txome")) {
    print(txome) 
} else {
    txome <- NULL
}

if (exists("run_camera")) {
    print(run_camera) 
} else {
    run_camera <- NULL
}


## Read metadata
msg0 <- try({
    metadata <- read.delim(metafile, header = TRUE, 
                           as.is = TRUE, sep = "\t");
    rownames(metadata) <- metadata$names;
    type <- metadata$type;
    if (length(unique(type)) == 1 | any(type %in% c("PE", "SE"))) {
        error("The type column in the metadata should have either PE or SE. \n")
    }
    }, 
    silent = TRUE)


msg1 <- try({
    if (unique(type) == "SE") {
        pt <- paste0(metadata$names, ".", fqsuffix, ".gz")
        
    } else {
        pt1 <- paste0(metadata$names, "_", fqext1, ".", fqsuffix, ".gz")
        pt2 <- paste0(metadata$names, "_", fqext2, ".", fqsuffix, ".gz")
        pt <- c(pt1, pt2)
    }
    
    lf <- file.path(fastqdir, pt)
    fe <- file.exists(lf)
    if (any(!fe)) {
        
        stop(fe[!fe], " are/is not available. \n")
    }
}, silent = TRUE)

msg2 <- try({
    fe <- file.exists(genome)
    if (!fe) {
        stop("The genonme file doesn't exist. \n")
    }
}, silent = TRUE)

msg3 <- try({
    fe <- file.exists(gtf)
    if (!fe) {
        stop("The gtf file doesn't exist. \n")
    }
}, silent = TRUE)

msg4 <- try({
    fe <- file.exists(txome)
    if (!fe) {
        stop("The txome file doesn't exist. \n")
    }
}, silent = TRUE)

# msg3 <- try({
#     if (run_camera == "True") {
#     library(msigdbr)
#         msigdbr_show_species()
#     }
# }, silent = TRUE)

msg6 <- try({
    if (exists("design")) {
        print(design)
    } else {
        design <- NULL
    }
}, silent = TRUE)


msg7 <- try({
    if (exists("contrast")) {
        contrast <- strsplit(gsub(" ","",contrast), ",")[[1]]
        print(contrast)
    } else {
        contrast <- NULL
    }
}, silent = TRUE)

## ---------------------------Test run -------------------------------


## Define design matrix
msg8 <- try({des <- model.matrix(as.formula(design), data = metadata)},
            silent = TRUE)

# Define contrasts
msg9 <- try({contrasts <- as.data.frame(makeContrasts(contrasts = contrast, 
                                                      levels = des))},
            silent = TRUE)

msgL <- list(msg0, msg1, msg2, msg3, msg4, 
             #msg5, 
             msg6, msg7, msg8, msg9)
isError <- lapply(msgL, FUN = function(x) {class(x) == "try-error"})
isError <- unlist(isError)

msg <- msgL[isError]

if (length(msg) > 0) {
    capture.output(writeLines(msg[[1]]), file = outFile)
} else {
    mylist <- list("Design matrix" = des, "Contrasts matrix" = contrasts)
    capture.output(mylist, file = outFile)
}



