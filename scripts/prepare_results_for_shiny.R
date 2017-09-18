args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(edgerres)
print(gtffile)
print(tx2gene)
print(metafile)
print(bigwigdir)
print(bwcolorvar)
print(outrds)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(rtracklayer))

options(ucscChromosomeNames = FALSE)

edgerres <- readRDS(edgerres)
tx2gene <- readRDS(tx2gene)

genes <- tx2gene %>% dplyr::select(-tx, -tx_biotype) %>% dplyr::distinct()

## -------------------------------------------------------------------------- ##
##                                 edgeR                                      ##
## -------------------------------------------------------------------------- ##
## edgeR result table with all contrasts ("wide")
edgeRwide <- edgerres$results
if (class(edgeRwide) == "list" && class(edgeRwide) != "data.frame") {
  for (i in seq_len(length(edgeRwide))) {
    edgeRwide[[i]]$F <- NULL
    colnames(edgeRwide[[i]]) <- paste0(colnames(edgeRwide[[i]]), ".", names(edgeRwide)[i])
    colnames(edgeRwide[[i]])[grep("logCPM", colnames(edgeRwide[[i]]))] <- "logCPM"
  }
  edgeRwide <- Reduce(function(...) dplyr::full_join(..., by = c("gene", "logCPM")), edgeRwide)
}

edgeRwide <- dplyr::left_join(edgeRwide, genes, by = "gene") %>%
  dplyr::select(gene, symbol, gene_biotype, logCPM, everything()) %>%
  dplyr::mutate(gene_biotype = factor(gene_biotype))

## edgeR result table for volcano plots ("long")
edgeRlong <- edgeRwide %>% tidyr::gather(typecontrast, value, -gene, -symbol, -gene_biotype, -logCPM) %>%
  dplyr::mutate(typecontrast = gsub("logFC.", "logFC_", typecontrast)) %>%
  dplyr::mutate(typecontrast = gsub("FDR.", "FDR_", typecontrast)) %>%
  dplyr::mutate(typecontrast = gsub("PValue.", "PValue_", typecontrast)) %>%
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
  genemodels <- subset(genemodels, type == "exon")
  
  genemodels
}

genemodels <- create_genemodels(gtffile)

## -------------------------------------------------------------------------- ##
##                          bigWig files + condition                          ##
## -------------------------------------------------------------------------- ##
## Vector with bigWig file names (relative to shiny app location) and condition information
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE)

bwfiles <- paste0("../", bigwigdir, "/", 
                  list.files(bigwigdir, pattern = "\\.bw$", full.names = FALSE))
names(bwfiles) <- gsub("\\.bw", "", basename(bwfiles))
condition <- sapply(names(bwfiles), function(w) {
  metadata[[bwcolorvar]][match(gsub("_trimmed_Aligned.sortedByCoord.out", "", w), metadata$ID)]
})

## -------------------------------------------------------------------------- ##
##                              edgeR - MDS                                   ##
## -------------------------------------------------------------------------- ##
## Data for PCA
dge <- edgerres$data
logcpms <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
mds <- limma::plotMDS(logcpms, top = 500, labels = NULL, pch = NULL, 
                      cex = 1, dim.plot = c(1, 2), ndim = 3, 
                      gene.selection = "common", 
                      xlab = NULL, ylab = NULL, plot = FALSE)$cmdscale.out
colnames(mds) <- c("MDS1", "MDS2", "MDS3")
mds <- as.data.frame(mds) %>% tibble::rownames_to_column(var = "ID") %>%
  dplyr::full_join(metadata)

## -------------------------------------------------------------------------- ##
##                                 Save                                       ##
## -------------------------------------------------------------------------- ##
saveRDS(list(results = list(edgeRwide = edgeRwide), 
             edgeRlong = edgeRlong, gene_models = genemodels,
             bw_files = bwfiles, condition = condition, pca = pca, 
             pcavar = pcavar, logcpms = logcpms, metadata = metadata, 
             gene_info = genes),
        file = outrds)


sessionInfo()
date()