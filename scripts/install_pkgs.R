args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(outtxt)

# =============== list all the packages =================
pkgs.use <- list(cran = c("dplyr", "ggplot2", "tidyr", "remotes"),
                 bioconductor = c("limma", "edgeR", "S4Vectors", "DRIMSeq", 
                                  "SingleCellExperiment", "tximeta"))

# install the packages if not installed
usePackage <- function(pkgs, defaultCRANmirror) {
    
    # install BiocManager package
    isBiocM <- "BiocManager" %in% installed.packages()[, 1]
    if (!isBiocM) {
        install.packages("BiocManager", repos = defaultCRANmirror)
        # install.packages("BiocManager")
    }
    
    # install the other packages
    pkgs <- unlist(pkgs.use, use.names = FALSE)
    isInstalled <- pkgs %in% installed.packages()[, 1]
    BiocManager::install(pkgs[!isInstalled],
                         update = FALSE, dependencies = TRUE)
    
    pkg.load <- lapply(pkgs, FUN = function(x) {
        x[!(x %in% installed.packages()[, "Package"])]
    })
    
    if (length(unlist(pkg.load)) == 0) {
        cat("All required packages are installed \n")
    } else {
        cat(unlist(pkg.load), ": failed to install")
    }
    
    suppressPackageStartupMessages(lapply(unlist(pkgs), library, character.only = TRUE))
    
    sink(outtxt)
    cat("packages loaded successfully: \n",
        unlist(pkgs.use)[unlist(pkgs) %in% loadedNamespaces()])
    sink()
}

paths <- .libPaths()
print(paths)

## install packages
usePackage(pkgs.use, defaultCRANmirror = "http://cran.at.r-project.org")
# usePackage(pkgs.use)
sessionInfo()
date()




