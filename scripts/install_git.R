args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(outtxt)


# =============== list all the packages =================
pkgs.use <- list(cran = c("dplyr","ggplot2", "tidyr", "devtools"),
                 bioconductor = c("limma", "edgeR", "ggtree",
                                  "S4Vectors", "DRIMSeq", 
                                  "SingleCellExperiment"),
                 github = c("tximeta"))

gitLink <- c("mikelove/tximeta")
names(gitLink) <- pkgs.use$github


# =================== load package =================
# install the packages if not installed
usePackage <- function(pkgs, gitLink, lib.personal,
                       defaultCRANmirror) {
    
    if (as.Date("2018-11-01") - Sys.Date()) {
    ## use this before 2018-10-31
    # --------------------------------------------------------------------------
    ## CRAN
    if(!is.null(pkgs[["cran"]])){
        ind.cran <- pkgs[["cran"]] %in% installed.packages()[, 1]
        if(sum(ind.cran)>0){
            cat("packages: ", pkgs[["cran"]][ind.cran], "exist in the environment \n")
        }
        if(!all(ind.cran)){
            install.packages(pkgs[["cran"]][!ind.cran], lib = lib.personal,
                             repos = defaultCRANmirror, quiet = TRUE)
            cat("packages: ", pkgs[["cran"]][!ind.cran], "are installing... \n")
        }
        
    }
    ## bioconductor
    if(!is.null(pkgs[["bioconductor"]])){
        ind.bioc <- pkgs[["bioconductor"]] %in% installed.packages()[, 1]
        if(sum(ind.bioc)>0){
            cat("packages: ", pkgs[["bioconductor"]][ind.bioc],
                "exist in the environment \n")
        }
        if(!all(ind.bioc)){
            cat("packages: ", pkgs[["bioconductor"]][!ind.bioc],
                "are installing... \n")
            source("https://bioconductor.org/biocLite.R")
            biocLite(pkgs[["bioconductor"]][!ind.bioc], suppressUpdates = TRUE,
                     lib = lib.personal)
        }
    }
    
    ## github (!!! remember to remove token)
    if(!is.null(pkgs[["github"]])){
        ind.git <- pkgs[["github"]] %in% installed.packages()[,1]
        if(sum(ind.git)>0){
            cat("packages: ", pkgs[["github"]][ind.git],
                "exist in the environment \n")
        }
        if(!all(ind.git)){
            options(unzip = "internal")
            for(i in seq_len(sum(!ind.git))){
                wgit <- which(!ind.git)
                cat("packages: ", pkgs[["github"]][wgit[i]],
                    "is installing \n")
                devtools::install_github(gitLink[pkgs[["github"]][wgit[i]]],
                                         lib = lib.personal,
                                         ref = "88271a3b57e7bafa93025fa09842d626abf036a1")
            }
        }
    }
    
    } else {
        
        ## use this after 2018-10-31
        # --------------------------------------------------------------------------
        pkgs <- unlist(pkgs.use, use.names = FALSE)
        isInstalled <- pkgs %in% installed.packages()[, 1]
        BiocManager::install(pkgs[!isInstalled], update = FALSE,
                             lib = lib.personal)
    }
    
    
    pkg.load <- lapply(pkgs, FUN = function(x){
        x[!(x %in% installed.packages()[,"Package"])]})
    
    if(length(unlist(pkg.load)) == 0){
        cat("All packages required are installed \n")
    }else{
        cat(unlist(pkg.load), ": failed to install")
    }
    
    # suppressPackageStartupMessages(lapply(unlist(pkgs), library, character.only = TRUE))
    
    sink(outtxt)
    cat("packages loaded successfully: \n",
        unlist(pkgs.use)[unlist(pkgs) %in% loadedNamespaces()])
    sink()  
    
}

paths <- .libPaths()
print(paths)

## install packages
usePackage(pkgs.use, gitLink = gitLink,
           lib.personal = .libPaths()[1],
           defaultCRANmirror = "http://cran.at.r-project.org")


sessionInfo()
date()




