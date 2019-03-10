args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
}

print(outtxt)

(mirror <- getOption("repos"))

## The function to install packages that are not installed
usePackage <- function(pkgs) {
    
    # install BiocManager package
    isBiocM <- "BiocManager" %in% installed.packages()[, 1]
    if (!isBiocM) {
        install.packages("BiocManager", repos = "http://cran.rstudio.com/",
                         lib = .libPaths()[1])
    }
    
    # install the other packages
    isInstalled <- pkgs %in% installed.packages()[, 1]
    BiocManager::install(pkgs[!isInstalled],
                         update = FALSE, dependencies = TRUE,
                         lib = .libPaths()[1])
    
    pkg.load <- lapply(pkgs, FUN = function(x) {
        x[!(x %in% installed.packages()[, "Package"])]
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

usePackage(pkgs = pkgs.use)


## Session info
sessionInfo()
date()

