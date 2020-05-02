args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

print(outtxt)
print(annotation)
print(organism)
print(ncores)

(mirror <- getOption("repos"))

## Function to install packages that are not installed
usePackage <- function(pkgs) {
    
    ## Install BiocManager package
    isBiocM <- "BiocManager" %in% installed.packages()[, 1]
    if (!isBiocM) {
        install.packages("BiocManager", repos = "http://cran.rstudio.com/",
                         lib = .libPaths()[1])
    }
    
    ## Check that Bioc is new enough
    if (BiocManager::version() < '3.11') {
      stop("Bioconductor release 3.11 or newer is required ", 
           "for this version of ARMOR.")
    }
    
    ## Install the other packages
    isInstalled <- pkgs %in% installed.packages(lib.loc = .libPaths()[1])[, 1]
    BiocManager::install(pkgs[!isInstalled],
                         update = FALSE, dependencies = TRUE,
                         lib = .libPaths()[1], Ncpus = as.integer(ncores))
    
    pkg.load <- lapply(pkgs, FUN = function(x) {
        x[!(x %in% installed.packages(.libPaths()[1])[, "Package"])]
    })
    
    if (length(unlist(pkg.load)) == 0) {
        cat("All required packages are installed \n")
    } else {
        cat(unlist(pkg.load), ": failed to install")
    }
    
    ## Test whether packages could be loaded successfully
    suppressPackageStartupMessages(
        lapply(pkgs, library, character.only = TRUE)
    )
    
    sink(outtxt)
    cat("packages loaded successfully: \n",
        pkgs[pkgs %in% loadedNamespaces()])
    sink()
}


paths <- .libPaths()
print(paths)

## Install packages
pkgs.use <- c("dplyr", "ggplot2", "tidyr", "remotes", "limma", "edgeR",
              "S4Vectors", "DRIMSeq", "SingleCellExperiment", "tximeta", "msigdbr")


if (annotation == "Gencode") {
  if (organism == "Homo_sapiens") {
    pkgs.extra = "org.Hs.eg.db"
  } else {
    pkgs.extra = "org.Mm.eg.db"
  }
  pkgs.use <- c(pkgs.use, pkgs.extra)
}
  

usePackage(pkgs = pkgs.use)


## Session info
sessionInfo()
date()

