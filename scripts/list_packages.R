args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## List all packages used in the analyses together with the version, by parsing
## the files in the "Routdir" directory. The results are written to the "outtxt"
## text file.

print(Routdir)
print(outtxt)

lf <- list.files(Routdir)
all_packages <- c()
for (f in lf) {
  x <- readLines(paste0(Routdir, "/", f))
  idx1 <- which(x == "other attached packages:")
  idx2 <- which(x == "loaded via a namespace (and not attached):")
  if (length(idx1) != 0 & length(idx2) != 0) {
    all_packages <- 
      unique(c(all_packages, 
               do.call(c, lapply((idx1 + 1):(idx2 - 2), function(i) {
                 grep("\\[", setdiff(setdiff(strsplit(x[i], " ")[[1]], " "), ""), 
                      value = TRUE, invert = TRUE)
               }))))
  }
}
write.table(sort(all_packages), file = outtxt, 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")