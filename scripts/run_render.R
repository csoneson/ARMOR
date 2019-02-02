args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## Mandatory arguments
print(se)
print(rmdtemplate)
print(outputdir)
print(outputfile)

## Arguments that are only used for some of the reports
if (exists("organism")) {
  print(organism)
} else {
  organism <- NULL
}

if (exists("gtffile")) {
  print(gtffile)
} else {
  gtffile <- NULL
}

if (exists("groupvar")) {
  print(groupvar)
} else {
  groupvar <- NULL
}

if (exists("bigwigdir")) {
  bigwigdir <- normalizePath(bigwigdir)
  print(bigwigdir)
} else {
  bigwigdir <- NULL
}

source("scripts/generate_report.R")

generateReport(se = se, organism = organism, gtffile = gtffile,
               bigwigdir = bigwigdir, groupvar = groupvar, 
               rmdTemplate = rmdtemplate, outputDir = outputdir, 
               outputFile = outputfile, forceOverwrite = TRUE,
               showCode = TRUE)