args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(edgerres)
print(gtffile)
print(tx2gene)
print(metafile)
print(bigwigdir)
print(groupvar)
print(outrds)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(reshape2))

options(ucscChromosomeNames = FALSE)

edgerres <- readRDS(edgerres)
tx2gene <- readRDS(tx2gene)

genes <- tx2gene %>% dplyr::select(-tx, -tx_biotype, -start, -end) %>% dplyr::distinct()

## -------------------------------------------------------------------------- ##
##                                 edgeR                                      ##
## -------------------------------------------------------------------------- ##
## edgeR result table with all contrasts ("wide")
edgeRwide <- edgerres$results
if (class(edgeRwide) == "list" && class(edgeRwide) != "data.frame") {
  for (i in seq_len(length(edgeRwide))) {
    idx <- which(colnames(edgeRwide[[i]]) %in% c("logFC", "F", "PValue", "FDR", "LR"))
    colnames(edgeRwide[[i]])[idx] <- paste0(colnames(edgeRwide[[i]])[idx], ".", names(edgeRwide)[i])
  }
  edgeRwide <- Reduce(function(...) dplyr::full_join(...), edgeRwide)
}

edgeRwide <- edgeRwide %>% 
  dplyr::select(gene, symbol, gene_biotype, logCPM, everything()) %>%
  dplyr::mutate(gene_biotype = factor(gene_biotype)) %>%
  dplyr::mutate(strand = factor(strand))

## edgeR result table for volcano plots ("long")
edgeRlong <- edgeRwide %>% 
  tidyr::gather(typecontrast, value, -gene, -symbol, -gene_biotype, -logCPM, 
                -chromosome, -strand) %>%
  dplyr::mutate(typecontrast = gsub("logFC\\.", "logFC_", typecontrast)) %>%
  dplyr::mutate(typecontrast = gsub("FDR\\.", "FDR_", typecontrast)) %>%
  dplyr::mutate(typecontrast = gsub("PValue\\.", "PValue_", typecontrast)) %>%
  dplyr::mutate(typecontrast = gsub("F\\.", "F_", typecontrast)) %>%
  dplyr::mutate(typecontrast = gsub("LR\\.", "LR_", typecontrast)) %>%
  tidyr::separate(typecontrast, into = c("dtype", "contrast"), sep = "_") %>%
  tidyr::spread(key = dtype, value = value) %>%
  dplyr::mutate(mlog10PValue = -log10(PValue))

## -------------------------------------------------------------------------- ##
##                             Gene models                                    ##
## -------------------------------------------------------------------------- ##
## Gene models from gtf
create_genemodels <- function(gtf_file) {
  genemodels <- rtracklayer::import(gtf_file)
  idx <- match(c("transcript_id", "gene_id", "exon_id"), colnames(mcols(genemodels)))
  colnames(mcols(genemodels))[idx] <- c("transcript", "gene", "exon")
  mcols(genemodels)$symbol <- mcols(genemodels)$transcript
  subset(genemodels, type == "exon")
}

if (!is.null(gtffile)) {
  genemodels <- create_genemodels(gtffile)
} else {
  genemodels <- NULL
}

## -------------------------------------------------------------------------- ##
##                          bigWig files + condition                          ##
## -------------------------------------------------------------------------- ##
## Vector with bigWig file names and condition information
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE)

if (!is.null(bigwigdir)) {
  bwfiles <- normalizePath(list.files(bigwigdir, pattern = "\\.bw$", full.names = TRUE))
  names(bwfiles) <- gsub("_Aligned.sortedByCoord.out.bw", "", basename(bwfiles))
  condition <- metadata[[groupvar]][match(names(bwfiles), metadata$ID)]
  names(condition) <- names(bwfiles)
  ordr <- order(condition)
  condition <- condition[ordr]
  bwfiles <- bwfiles[ordr]
} else {
  bwfiles <- condition <- NULL
}

## -------------------------------------------------------------------------- ##
##                              edgeR - MDS                                   ##
## -------------------------------------------------------------------------- ##
## Data for MDS
dge <- edgerres$data
logcpms <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
mds <- limma::plotMDS(logcpms, top = 500, labels = NULL, pch = NULL, 
                      cex = 1, dim.plot = c(1, 2), ndim = min(7, ncol(logcpms) - 1), 
                      gene.selection = "common", 
                      xlab = NULL, ylab = NULL, plot = FALSE)$cmdscale.out
colnames(mds) <- paste0("MDS", 1:min(7, ncol(logcpms) - 1))
mds <- as.data.frame(mds) %>% tibble::rownames_to_column(var = "ID") %>%
  dplyr::full_join(metadata)

logcpms <- reshape2::melt(as.matrix(logcpms)) %>%
  dplyr::rename(gene = Var1, sample = Var2) %>%
  dplyr::mutate(group = metadata[match(sample, metadata$ID), groupvar])

## -------------------------------------------------------------------------- ##
##                                 Save                                       ##
## -------------------------------------------------------------------------- ##
saveRDS(list(wideResults = list(edgeR = edgeRwide), 
             longResults = list(edgeR = edgeRlong), 
             geneModels = genemodels,
             bwFiles = bwfiles, 
             bwCond = condition, 
             dimReds = list(MDS = mds),
             abundances = list(logCPM = logcpms), 
             geneInfo = genes),
        file = outrds)


sessionInfo()
date()