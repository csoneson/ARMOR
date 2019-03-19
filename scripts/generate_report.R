suppressPackageStartupMessages({
    library(rmarkdown)
})

#' Generate report
#'
#' Generate a report based on a Rmarkdown template file.
#' 
#' @param se,gtffile,bigwigdir,genesets,organism,design,contrast Arguments that
#'   are passed to the provided Rmarkdown template
#' @param rmdTemplate Path to a .Rmd template file.
#' @param outputFile File name of the output report. The file name extension
#'   must be either \code{.html} or \code{.pdf}, and consistent with the value
#'   of \code{outputFormat}.
#' @param outputDir Path to the output directory where the report will be
#'   generated.
#' @param outputFormat The format of the output report. Either
#'   \code{"html_document"} or \code{"pdf_document"}. The file name extension of
#'   \code{outputFile} must be consistent with this choice.
#' @param showCode Logical, whether to display the R code in the report.
#' @param forceOverwrite Logical, whether to force overwrite an existing report
#'   with the same name in the output directory.
#' @param knitrProgress Logical, whether to display the progress of \code{knitr}
#'   when generating the report.
#' @param quiet Logical, whether to show progress messages.
#' @param ignorePandoc Logical, determines what to do if \code{pandoc} or
#'   \code{pandoc-citeproc} is missing (if \code{Sys.which("pandoc")} or
#'   \code{Sys.which("pandoc-citeproc")} returns ""). If \code{ignorePandoc} is
#'   TRUE, only a warning is given. The figures will be generated, but not the
#'   final report. If \code{ignorePandoc} is FALSE (default), the execution
#'   stops immediately.
#' @param ... Other arguments that will be passed to \code{rmarkdown::render}.
#'
#' @author Charlotte Soneson
#'
#' @details When the function is called, an .Rmd template file will be copied
#'   into the output directory, and \code{rmarkdown::render} will be called to
#'   generate the final report. If there is already a .Rmd file with the same
#'   name in the output directory, the function will raise an error and stop, to
#'   avoid overwriting the existing file. The reason for this behaviour is that
#'   the copied template in the output directory will be deleted once the report
#'   is generated.
#'
#' @export
#'
#' @importFrom rmarkdown render
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom methods is
#' @import dplyr
#'
#' @return Generates a summary report in the \code{outputDir} directory, and
#'   returns (invisibly) the name of the generated report.
#'

generateReport <- function(se, gtffile = NULL, organism = NULL, 
                           bigwigdir = NULL, design = NULL, genesets = NULL, 
                           contrast = NULL, rmdTemplate, outputFile, 
                           ncores = NULL,
                           outputDir = "./", outputFormat = NULL, 
                           showCode = FALSE, forceOverwrite = FALSE, 
                           knitrProgress = FALSE, quiet = FALSE, 
                           ignorePandoc = FALSE, ...) {
    ## This function was inspired by code from Nicholas Hamilton, provided at
    ## http://stackoverflow.com/questions/37097535/generate-report-in-r
    
    ## If possible, set output format based on the extension of outputFile, if
    ## the output format is not provided
    if (is.null(outputFormat)) {
        if (tools::file_ext(outputFile) == "pdf") {
            outputFormat <- "pdf_document"
        } else {
            outputFormat <- "html_document"
        }
    }
    
    ## Check if pandoc and pandoc-citeproc are available
    if (Sys.which("pandoc") == "") {
        if (ignorePandoc) {
            ## If ignorePandoc is TRUE, just give a warning
            warning("pandoc is not available! ",
                    "The final report will not be generated.")
        } else {
            ## If ignorePandoc is FALSE, stop
            stop("pandoc is not available!")
        }
    }
    if (Sys.which("pandoc-citeproc") == "") {
        if (ignorePandoc) {
            ## If ignorePandoc is TRUE, just give a warning
            warning("pandoc-citeproc is not available! ",
                    "The final report will not be generated.")
        } else {
            ## If ignorePandoc is FALSE, stop
            stop("pandoc-citeproc is not available!")
        }
    }
    
    ## ---------------------------------------------------------------------- ##
    ## --------------------- Check input arguments -------------------------- ##
    ## ---------------------------------------------------------------------- ##
    
    ## ------------------------ outputFormat -------------------------------- ##
    ## Raise an error if outputFormat is not one of the allowed
    if (!(outputFormat %in% c("pdf_document", "html_document"))) {
        stop("The provided outputFormat is currently not supported. Please ",
             "use either 'html_document' or 'pdf_document'.", call. = FALSE)
    }
    
    ## Raise an error if the output format and file name extension don't match
    if (outputFormat != paste0(tools::file_ext(outputFile), "_document")) {
        stop(paste0("File name extension of outputFile doesn't agree with the ",
                    "outputFormat, should be .",
                    gsub("_document$", "", outputFormat)), call. = FALSE)
    }
    
    ## ----------------------- input directory ------------------------------ ##
    ## se must be a character string of length 1, and point to an existing rds file
    if (!is(se, "character") || length(se) != 1) {
        stop("se must be a character string")
    }
    if (!file.exists(se)) {
        stop("The indicated se object does not exist")
    }
    se <- readRDS(se)
    
    ## organism
    if (!is.null(organism)) {
        if (!is(organism, "character") || length(organism) != 1) {
            stop("organism must be a character string")
        }
        organism <- gsub("_", " ", organism)
        if (!is.null(genesets) && 
            !organism %in% msigdbr::msigdbr_show_species()) {
            stop("organism must be one of the organisms listed in ",
                 "msigdbr::msigdbr_show_species()")
        }
    }
    
    ## design
    if (!is.null(design)) {
        if (!is(design, "character") || length(design) != 1) {
            stop("design must be a character string")
        }
    }
    
    ## contrasts
    if (!is.null(contrast)) {
        if (!is(contrast, "character")) {
            stop("contrast must be a character string")
        }
    }

    ## ncores
    if (!is.null(ncores)) {
        if (!is(ncores, "numeric")) {
            stop("ncores must be numeric")
        }
    }
    
    ## genesets
    if (!is.null(genesets)) {
        if (!is(genesets, "character")) {
            stop("genesets must be a character string")
        }
    }
    
    ## gtffile 
    if (!is.null(gtffile)) {
        if (!is(gtffile, "character") || length(gtffile) != 1) {
            stop("gtffile must be a character string")
        }
        if (!file.exists(gtffile)) {
            stop("The indicated gtffile does not exist")
        }
        genemodels <- rtracklayer::import(gtffile)
    }
    
    ## bigwigdir
    if (!is.null(bigwigdir)) {
        if (!is(bigwigdir, "character") || length(bigwigdir) != 1) {
            stop("bigwigdir must be a character string")
        }
        if (!file.exists(bigwigdir)) {
            stop("The indicated bigwigdir does not exist")
        }
    }
    
    ## ------------------------- output files ------------------------------- ##
    outputReport <- file.path(outputDir, basename(outputFile))
    outputRmd <- file.path(
        outputDir,
        paste0(tools::file_path_sans_ext(basename(outputFile)), ".Rmd"))
    
    ## Report
    if (file.exists(outputReport)) {
        if (!forceOverwrite) {
            stop("The file ", outputReport,
                 " already exists. Please remove or rename the file, provide ",
                 "another value of outputFile, or set forceOverwrite = TRUE.",
                 call. = FALSE)
        } else {
            if (!quiet) {
                warning("The file ", outputReport,
                        " already exists and will be overwritten, since ",
                        "forceOverwrite = TRUE.", immediate. = TRUE,
                        call. = FALSE)
            }
        }
    }
    
    ## ------------------------- Rmd template ------------------------------- ##
    ## Path to the template file
    templateFile <- rmdTemplate
    if (file.exists(templateFile)) {
        if (file.exists(outputRmd)) {
            if (!forceOverwrite) {
                stop("There is already an .Rmd file ", outputRmd,
                     ". Please remove or rename this file, or choose another ",
                     "outputFile name.", call. = FALSE)
            } else {
                warning("There is already an .Rmd file ", outputRmd, 
                        ". That file will be renamed with a suffix '_conflicting'",
                        ", a time stamp and a random sequence. If you did not ", 
                        "explicitly create this file, it can be removed.", 
                        call. = FALSE)
                file.rename(from = outputRmd, 
                            to = paste0(outputRmd, "_conflicting_", Sys.Date(), "_", 
                                        round(1e6*runif(1))))
            }
        }
        file.copy(from = templateFile, to = outputRmd, overwrite = FALSE)
    } else {
        stop("The Rmd template file ", templateFile, " does not exist.",
             call. = FALSE)
    }
    
    ## ---------------------------------------------------------------------- ##
    ## ----------------------- Process the arguments ------------------------ ##
    ## ---------------------------------------------------------------------- ##
    
    args <- list(...)
    args$input <- outputRmd
    args$output_format <- outputFormat
    args$output_file <- outputFile
    args$quiet <- !knitrProgress
    
    ## ---------------------------------------------------------------------- ##
    ## ------------------------ Render the report --------------------------- ##
    ## ---------------------------------------------------------------------- ##
    
    outputFile <- do.call("render", args = args)
    
    ## ---------------------------------------------------------------------- ##
    ## --------------------- Remove temporary file -------------------------- ##
    ## ---------------------------------------------------------------------- ##
    
    file.remove(outputRmd)
    
    invisible(outputFile)
}
