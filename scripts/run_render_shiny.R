args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
}

print(se)
print(gtffile)
print(groupvar)
print(rmdtemplate)
print(outputdir)
print(outputfile)
if (exists("bigwigdir")) print(bigwigdir)

source("scripts/generate_report.R")

generateReport(se = se, gtffile = gtffile,
               bigwigdir = if (exists("bigwigdir")) normalizePath(bigwigdir) else NULL, 
               groupvar = groupvar, 
               rmdTemplate = rmdtemplate, outputDir = outputdir, 
               outputFile = outputfile, forceOverwrite = TRUE,
               showCode = TRUE)
