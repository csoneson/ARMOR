args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(se)
print(organism)
print(rmdtemplate)
print(outputdir)
print(outputfile)

source("scripts/generate_report.R")

generateReport(se = se, organism = organism, 
               rmdTemplate = rmdtemplate, outputDir = outputdir, 
               outputFile = outputfile, forceOverwrite = TRUE,
               showCode = TRUE)
