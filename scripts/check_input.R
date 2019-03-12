
args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}


## ----------------The input arguments------------------------------------------
if (exists("metafile")) {
    print(metafile) 
    } else {
        metafile <- NULL
    }

if (exists("organism")) {
    print(organism)
} else {
    organism <- NULL
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
    if(!file.exists(metafile)) {
        error("The metafile, ", metafile, ", does not exist.\n")
    } else {
        metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t");
        if(!all(c("names","type") %in% colnames(metadata)))
            stop(paste0("ERROR: 'names' and 'type' columns must exist in ", metafile))
        rownames(metadata) <- metadata$names;
        utype <- unique(metadata$type);
        if (length(utype) == 1 & any(utype %in% c("PE", "SE"))) {
            type <- metadata$type
        } else{
            stop("ERROR: 'type' column in the metadata file must be PE or SE. \n")
        }
    }
}, silent = TRUE)


msg1 <- try({
    if (utype == "SE") {
        pt <- paste0(metadata$names, ".", fqsuffix, ".gz")
    } else {
        pt1 <- paste0(metadata$names, "_", fqext1, ".", fqsuffix, ".gz")
        pt2 <- paste0(metadata$names, "_", fqext2, ".", fqsuffix, ".gz")
        pt <- c(pt1, pt2)
    }
    lf <- file.path(fastqdir, pt)
    fe <- file.exists(lf)
    if (any(!fe)) {
        stop(paste0("ERROR: ", paste(lf[!fe], collapse=" "), " are/is not available.\n"))
    }
}, silent = TRUE)

print(lf)
print(fe)

msg2 <- try({
    fe <- file.exists(genome)
    if (!fe) {
        stop(paste0("ERROR: The 'genome' file, ", genome, ", doesn't exist.\n"))
    }
}, silent = TRUE)

msg3 <- try({
    fe <- file.exists(gtf)
    if (!fe) {
        stop(paste0("ERROR: The 'gtf' file, ", gtf, ", doesn't exist.\n"))
    }
}, silent = TRUE)

msg4 <- try({
    fe <- file.exists(txome)
    if (!fe) {
        stop(paste0("ERROR: The 'txome' file, ", txome, ", doesn't exist.\n"))
    }
}, silent = TRUE)

msg5 <- try({
    if (run_camera == "True")
      if (require("msigdbr")) {
          if (!(gsub("_"," ",organism) %in% msigdbr::msigdbr_show_species()))
              stop(paste0("ERROR: '", gsub("_"," ",organism), "' not found in 'msigdbr::msigdbr_show_species()' database; fix the organism or set 'run_camera: False'"))
      } else {
          stop("Cannot check 'organism': msigdbr package not available; run 'snakemake [--use-conda] setup' before 'snakemake [--use-conda] checkinputs'")
      }
    }, silent = TRUE)

msg6 <- try({
    if (exists("design")) {
        print(design)
    } else {
        stop("ERROR: no 'design' specified; please specify one in the config file")
    }
}, silent = TRUE)


msg7 <- try({
    if (exists("contrast")) {
        contrast <- strsplit(gsub(" ","",contrast), ",")[[1]]
        print(contrast)
    } else {
        stop("ERROR: no 'contrast' specified; please specify one in the config file")
    }
}, silent = TRUE)

## Define design matrix
msg8 <- try({
    des <- model.matrix(as.formula(design), data = metadata)
}, silent = TRUE)
if(is(msg8, "try-error"))
    msg8 <- try({
        stop("ERROR in 'design' value: ", design)
    }, silent=TRUE)


# Define contrasts
msg9 <- try({
    have_edgeR <<- FALSE
    if (require("edgeR")) {
        have_edgeR <<- TRUE
        contrasts <- as.data.frame(makeContrasts(contrasts = contrast, 
                                                 levels = des))
    } else {
        stop("Cannot check 'contrast', since the edgeR package is not available; run 'snakemake [--use-conda] setup' before 'snakemake [--use-conda] checkinputs'")
    }
}, silent = TRUE)
if(is(msg9, "try-error") && have_edgeR)
    msg9 <- try({
        stop("ERROR in specified 'contrast' (n.b., could be due to 'design'): ", paste0(contrast, collapse=","))
    }, silent=TRUE)

msgL <- list(msg0, msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9)
isError <- sapply(msgL, FUN = function(x) {class(x) == "try-error"})
msg <- msgL[isError]
print(msg)

if (length(msg) > 0) {
    for(i in seq_len(length(msg))) {
      m <- trimws(gsub("Error in try({ :", "", msg[[i]], fixed=TRUE))
      capture.output(writeLines(m), file = outFile, append = !(i==1))
    }
    stars <- paste(strrep("*", 84), "\n", strrep("*", 84), sep="")
    xmsg <- paste("check for the error message above and fix the config.yaml or one of it's components.", sep="")
    capture.output(writeLines(stars), file = outFile, append = TRUE)
    capture.output(writeLines(xmsg), file = outFile, append = TRUE)
    capture.output(writeLines(stars), file = outFile, append = TRUE)
} else {
    mylist <- list("Design matrix" = des, "Contrasts matrix" = contrasts)
    capture.output(mylist, file = outFile)
    stars <- paste(strrep("*", 19), "\n", strrep("*", 19), sep="")
    xmsg <- paste("No errors detected.", sep="")
    capture.output(writeLines(stars), file = outFile, append = TRUE)
    capture.output(writeLines(xmsg), file = outFile, append = TRUE)
    capture.output(writeLines(stars), file = outFile, append = TRUE)
}



