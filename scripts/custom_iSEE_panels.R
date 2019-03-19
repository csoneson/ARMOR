suppressPackageStartupMessages({
    library(ggplot2)
    library(rtracklayer)
    library(iSEE)
})

options(ucscChromosomeNames = FALSE)
prepareGtf <- function(gtf) {
    gtf <- rtracklayer::import(gtf)
    
    ## Set appropriate column names
    idx <- match(c("transcript_id", "gene_id", "exon_id"),
                 colnames(S4Vectors::mcols(gtf)))
    colnames(S4Vectors::mcols(gtf))[idx] <- c("transcript", "gene", "exon")
    if (!("gene_name" %in% colnames(S4Vectors::mcols(gtf)))) {
        gtf$gene_name <- gtf$gene
    }
    
    ## Keep only exons
    gtf <- BiocGenerics::subset(gtf, type == "exon")
    
    ## Strip version numbers from gene and transcript IDs if they exist
    gtf$transcript <- gsub("\\.[0-9]+$", "", gtf$transcript)
    gtf$gene <- gsub("\\.[0-9]+$", "", gtf$gene)
    
    gtf
}

customGviz <- function(se, rows, columns, bigwig_files="", bigwig_names="", 
                       bigwig_condition="", granges="",
                       chr="", start="", end="", showgene="") {
    options(ucscChromosomeNames = FALSE)
    
    ## ---------------------------------------------------------------------- ##
    ## Pre-flight checks
    ## ---------------------------------------------------------------------- ##
    ## Must have at least one of bigwig_files and granges
    if (bigwig_files == "" && granges == "") {
        return(NULL)
    }
    
    ## If no names are given, assign names to bigwig files
    if (bigwig_files != "" && bigwig_names == "") {
        bigwig_names <- paste(paste0("S", seq_along(strsplit(bigwig_files, ",")[[1]])),
                              collapse = ",")
    }
    
    ## If granges file does not exist, don't show annotation
    if (!file.exists(granges)) {
        granges <- ""
    }
    
    ## If granges file does not exist, the viewing region must be set
    if (granges == "" && (chr == "" || start == "" || end == "")) {
        return(NULL)
    }
    
    ## Convert start and end positions to numeric values
    if (start != "") {
        start <- as.numeric(start)
    }
    if (end != "") {
        end <- as.numeric(end)
    }
    
    ## ---------------------------------------------------------------------- ##
    ## Prepare the annotation
    ## ---------------------------------------------------------------------- ##
    if (granges != "") {
        ## Read the GRanges object
        if (caching$granges == granges && !is.null(caching$gr0)) {
            gr0 <- caching$gr0
        } else {
            caching$gr0 <- readRDS(granges)
            caching$granges <- granges
            gr0 <- caching$gr0
        }
        
        ## Subset the GRanges object depending on the input
        ## If rows has length 1, overwrite any provided showgene
        if (length(rows) == 1) {
            showgene <- rows
        } 
        
        ## Strip version number from the gene of interest if it exists
        showgene <- gsub("\\.[0-9]+$", "", showgene)
        
        if (showgene == "" && (chr == "" || is.na(start) || is.na(end))) {
            return(NULL)
        }
        
        ## If a gene has been defined (either via rows or via showgene), set the 
        ## viewing range accordingly
        if (showgene != "") {
            gr <- BiocGenerics::subset(gr0, tolower(gene) == tolower(showgene) | 
                                           tolower(gene_name) == tolower(showgene))
            ## Select only one gene if there are many with the same name
            gr <- BiocGenerics::subset(gr, gene == gene[1])
            chr <- unique(GenomeInfoDb::seqnames(gr))
            start <- min(BiocGenerics::start(gr))
            end <- max(BiocGenerics::end(gr))
        } else {
            gr <- gr0[IRanges::overlapsAny(
                gr0,
                GenomicRanges::GRanges(seqnames = chr,
                                       ranges = IRanges::IRanges(start = start,
                                                                 end = end),
                                       strand = "*")), ]
        }
        
        ## Other features in the region
        gro <- gr0[IRanges::overlapsAny(
            gr0,
            GenomicRanges::GRanges(seqnames = chr,
                                   ranges = IRanges::IRanges(start = start,
                                                             end = end),
                                   strand = "*"))]
        gro <- gro[!(S4Vectors::`%in%`(gro, gr))]
        
        grtr <- Gviz::GeneRegionTrack(gr, showId = TRUE, col = NULL, fill = "gray80",
                                      name = "Genes", col.title = "black")
        grtr2 <- Gviz::GeneRegionTrack(gro, showId = TRUE, col = "black", fill = "white",
                                       name = "", col.title = "black")
    } else {
        gr <- gro <- grtr <- grtr2 <- NULL
    }    
    
    ## ---------------------------------------------------------------------- ##
    ## Set title and viewing range
    ## ---------------------------------------------------------------------- ##
    ## Define the title for the plot
    if (showgene != "" && !is.null(gr)) {
        if (all(gr$gene == gr$gene_name)) {
            plot_title <- unique(gr$gene)
        } else {
            plot_title <- unique(paste0(gr$gene, " (", gr$gene_name, ")"))
        }
    } else {
        plot_title <- paste0(chr, ":", start, "-", end)
    }
    
    ## Set min and max coord for the plot (add some padding to each side)
    minCoord <- start - 0.15*(end - start)
    maxCoord <- end + 0.05*(end - start)
    
    ## ---------------------------------------------------------------------- ##
    ## Prepare bigWig files
    ## ---------------------------------------------------------------------- ##
    ## Reformat bigWig file paths and names (provided to the function as 
    ## character strings)
    if (bigwig_files != "") {
        bigwig_files <- strsplit(bigwig_files, ",")[[1]]
        bigwig_names <- strsplit(bigwig_names, ",")[[1]]
        if (bigwig_condition != "") {
            bigwig_condition <- strsplit(bigwig_condition, ",")[[1]]
            names(bigwig_condition) <- bigwig_names
        }
        names(bigwig_files) <- bigwig_names
        
        ## ---------------------------------------------------------------------- ##
        ## Define colors if bigwig_condition is provided
        ## ---------------------------------------------------------------------- ##
        ## Define colors for coverage tracks
        color_list <- rep(c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D", "#F7EE55",
                            "#90C987", "#777777", "#E8601C", "#1965B0", "#882E72",
                            "#F6C141", "#4EB265", "#CAEDAB"), 
                          ceiling(length(unique(bigwig_condition))/13))
        
        if (length(bigwig_condition) > 1 || bigwig_condition != "") {
            usecol <- color_list[match(bigwig_condition, 
                                       unique(bigwig_condition))]
        } else {
            usecol <- rep("gray", length(bigwig_files))
        }
        names(usecol) <- bigwig_names
        
        ## ------------------------------------------------------------------ ##
        ## Show only selected sample(s)
        ## ------------------------------------------------------------------ ##
        ## If columns is specified, subset bigwig files
        if (!is.null(columns)) {
            bigwig_files <- bigwig_files[columns]
            bigwig_condition <- bigwig_condition[columns]
            usecol <- usecol[columns]
        }
        
        ## ------------------------------------------------------------------ ##
        ## Prepare final plot
        ## ------------------------------------------------------------------ ##
        ## Set up coverage tracks
        tracks <- lapply(seq_along(bigwig_files), function(i) {
            assign(paste0("covtr", i), 
                   Gviz::DataTrack(range = bigwig_files[i],
                                   type = "histogram",
                                   name = names(bigwig_files)[i],
                                   col.title = "black",
                                   fill = usecol[i],
                                   col = usecol[i],
                                   col.histogram = usecol[i],
                                   fill.histogram = usecol[i]))
        })
    } else {
        tracks <- NULL
    }
    
    ## Add genome axis track
    tracks <- c(tracks, Gviz::GenomeAxisTrack(), grtr, grtr2)
    
    ## Plot tracks
    Gviz::plotTracks(tracks, chromosome = chr, from = minCoord, 
                     to = maxCoord, main = plot_title, 
                     transcriptAnnotation = "transcript",
                     min.width = 0, min.distance = 0, collapse = FALSE)
}

customVolcano <- function(se, rows, columns, contrasts) {
    contrasts <- strsplit(contrasts, ",")[[1]]
    tmp <- do.call(plyr::rbind.fill, lapply(contrasts, function(w) {
        x <- data.frame(rowData(se)[, grep(paste0("^", w, ":"), 
                                           colnames(rowData(se))), 
                                    drop = FALSE], check.names = FALSE)
        colnames(x) <- gsub(paste0("^", w, ":"), "", colnames(x))
        x$contrast <- w
        x$feature <- rownames(x)
        x
    }))
    ggplot(tmp, aes(x = logFC, y = mlog10PValue)) + 
        geom_point(alpha = 0.3) + facet_grid(~ contrast) + 
        theme_bw() + ylab("-log10(PValue)")
}

# Set up a cache for the GRanges object
caching <- new.env()

# gtf <- prepareGtf("example_data/reference/Homo_sapiens.GRCh38.93.1.1.10M.gtf")
# saveRDS(gtf, file = "example_data/reference/Homo_sapiens.GRCh38.93.1.1.10M.granges.rds")
# 
# cdp <- customDataPlotDefaults(sce, 2)
# cdp$Function <- c("customGviz", "customVolcano")
# cdp$Arguments <- c("bigwig_files example_data/output/STARbigwig/SRR1039508_Aligned.sortedByCoord.out.bw,example_data/output/STARbigwig/SRR1039509_Aligned.sortedByCoord.out.bw,example_data/output/STARbigwig/SRR1039512_Aligned.sortedByCoord.out.bw,example_data/output/STARbigwig/SRR1039513_Aligned.sortedByCoord.out.bw\nbigwig_names SRR1039508,SRR1039509,SRR1039512,SRR1039513\nbigwig_condition Untreated,Dexamethasone,Untreated,Dexamethasone\ngranges example_data/reference/Homo_sapiens.GRCh38.93.1.1.10M.granges.rds\nchr 1\nstart 6.1e6\nend 6.2e6\nshowgene DDX11L1", 
#                    "contrasts cellineN61311-cellineN052611")
# 
# iSEE(sce, 
#      customDataArgs = cdp, 
#      customDataFun = list(customGviz = customGviz, 
#                           customVolcano = customVolcano))

