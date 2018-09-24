args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
}

print(edgerres)
print(gtffile)
print(metafile)
print(bigwigdir)
print(groupvar)
print(outList)
print(outSCE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(limma)
  library(edgeR)
  library(rtracklayer)
  library(reshape2)
  library(SingleCellExperiment)
  library(S4Vectors)
})

options(ucscChromosomeNames = FALSE)
edgerres <- readRDS(edgerres)

genes <- edgerres$data$genes %>%
    dplyr::select(gene = gene_id, symbol, gene_biotype, chromosome = seqnames, strand) %>%
    dplyr::mutate(strand = ifelse( strand == "-", -1,1)) %>%
    dplyr::distinct()

## -------------------------------------------------------------------------- ##
##                                 edgeR                                      ##
## -------------------------------------------------------------------------- ##
## edgeR result table with all contrasts ("wide")
edgeRwide <- edgerres$results
if (is(edgeRwide, "list")) {
    for (i in seq_len(length(edgeRwide))) {
        idx <- which(colnames(edgeRwide[[i]]) %in% c("logFC", "F", "PValue", "FDR", "LR"))
        colnames(edgeRwide[[i]])[idx] <- 
          paste0(colnames(edgeRwide[[i]])[idx], ".", names(edgeRwide)[i])
    }
    edgeRwide <- Reduce(function(...) dplyr::full_join(...), edgeRwide)
}


edgeRwide <- edgeRwide %>%
    dplyr::select(-start, -end, -width, -gene_name, -entrezid, -seq_coord_system) %>%
    dplyr::select(gene = gene_id, symbol, gene_biotype, logCPM, 
                  chromosome = seqnames, strand, everything()) %>%
    dplyr::mutate(strand = factor(ifelse(strand == "-", -1,1))) %>%
    dplyr::mutate(gene_biotype = factor(gene_biotype))


## edgeR result table for volcano plots ("long")
selC <- setdiff(colnames(edgeRwide), c("gene", "symbol", "gene_biotype",
                                       "logCPM", "chromosome", "strand"))
sepSel <- c(".", "_",":", "@", "%", "&")
chs <- sepSel[!unlist(lapply(seq_along(sepSel),
                             FUN = function(x){
                                 any(grepl(sepSel[x], selC))}))][1]

edgeRlong <- edgeRwide %>%
    tidyr::gather(typecontrast, value, -gene, -symbol, -gene_biotype, -logCPM,
                  -chromosome, -strand) %>%
    dplyr::mutate(typecontrast = gsub("logFC\\.", paste("logFC",chs,sep=""),
                                      typecontrast)) %>%
    dplyr::mutate(typecontrast = gsub("FDR\\.",  paste("FDR",chs,sep=""),
                                      typecontrast)) %>%
    dplyr::mutate(typecontrast = gsub("PValue\\.", paste("PValue",chs,sep=""), typecontrast)) %>%
    dplyr::mutate(typecontrast = gsub("F\\.", paste("F",chs,sep=""), typecontrast)) %>%
    dplyr::mutate(typecontrast = gsub("LR\\.", paste("LR",chs,sep=""), typecontrast)) %>%
    tidyr::separate(typecontrast, into = c("dtype", "contrast"), sep = chs) %>%
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
    condition <- metadata[[groupvar]][match(names(bwfiles), metadata$names)]
    names(condition) <- names(bwfiles)
    ordr <- order(condition)
    condition <- condition[ordr]
    bwfiles <- bwfiles[ordr]
} else {
    # bwfiles <- condition <- NULL
    bwfiles <- condition <- NA
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
mds <- as.data.frame(mds) %>% tibble::rownames_to_column(var = "names") %>%
    dplyr::full_join(metadata, by = "names")

logcpms <- reshape2::melt(as.matrix(logcpms)) %>%
    dplyr::rename(gene = Var1, sample = Var2) %>%
    dplyr::mutate(group = metadata[match(sample, metadata$names), groupvar])

## -------------------------------------------------------------------------- ##
##                             export - SCE  & SE                             ##
## -------------------------------------------------------------------------- ##
## row data
# information about genes
rData <- edgerres$data$genes %>%
    dplyr::distinct() %>%
    dplyr::arrange(match(gene_id, rownames(dge)))
rData <- S4Vectors::DataFrame(rData)

# output from edgeR
typeContrast <- unique(edgeRlong$contrast)
resList <- lapply(seq_along(typeContrast), FUN = function(x) {
    # rows kept in the analysis
    data.x <- edgeRlong %>%
        dplyr::filter(contrast == typeContrast[x]) %>%
        dplyr::select(-symbol, -gene_biotype, -logCPM,
                      -chromosome, -strand, -contrast)
    # add rows for those filtered out
    geneOUT <- setdiff(rownames(dge), data.x$gene)
    if (length(geneOUT) > 0) {
        data.x <- data.frame(matrix(NA, nrow = length(geneOUT),
                                    ncol = ncol(data.x),
                                    dimnames = list(geneOUT, colnames(data.x))
        )) %>%
            dplyr::mutate(gene = geneOUT) %>%
            dplyr::bind_rows(data.x) %>%
            dplyr::arrange(match(gene, rownames(dge)))
        
    } else {
        data.x <- data.x
    }
    
    # use gene column as rownames
    rownames(data.x) <- data.x$gene
    data.f <- data.x %>%
        dplyr::select(-gene)
    data.f <- S4Vectors::DataFrame(data.f)
    
    return(data.f)
})
names(resList) <- typeContrast

# put edgeR output in the rowData
for (i in seq_along(resList)) {
    tc.i <- typeContrast[i]
    rData[[tc.i]] <- resList[[i]]
}

## column data
# bwFiles and bwCond
cData1 <- data.frame(names = colnames(dge),
                     bwFiles = bwfiles[colnames(dge)],
                     bwCond = condition[colnames(dge)]) %>%
    dplyr::arrange(match(names, colnames(dge)))
# metadata
cData2 <- metadata %>%
    dplyr::arrange(match(names, colnames(dge)))
cData <- dplyr::bind_cols(cData1, cData2)
cData <- DataFrame(cData)

## assays data
logCPM_count <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
aData <- list(rawCount = dge$counts,
              logCPM = logCPM_count[rownames(dge), ])

## low dimensional representations
reducedData <- mds %>%
    dplyr::arrange(match(names, colnames(dge))) %>%
    dplyr::select(-one_of(colnames(metadata)))

reducedData <- as.matrix(reducedData)

sce <- SingleCellExperiment(assays = aData, rowData = rData,
                            colData = cData,
                            metadata = list(geneModels = genemodels,
                                            geneInfo = genes),
                            reducedDims = SimpleList(dimReds = reducedData))

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
        file = outList)

saveRDS(sce,
        file = outSCE)

sessionInfo()
date()
