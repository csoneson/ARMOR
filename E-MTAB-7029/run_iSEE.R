sce <- readRDS("E-MTAB-7029/output/outputR/shiny_sce.rds")
sce <- sce$sce_gene
rownames(sce) <- paste0(rowData(sce)$gene_id, "__", rowData(sce)$symbol)
rowData(sce) <- tidyr::unnest(as.data.frame(rowData(sce)))
rowData(sce)$edgeR.conditiond4Tcf__chir.conditiond4Tcf__unstim.mlog10PValue <- 
    -log10(rowData(sce)$edgeR.conditiond4Tcf__chir.conditiond4Tcf__unstim.PValue)

reddim <- redDimPlotDefaults(sce, 5)
reddim$PointSize <- 5
reddim$ColorBy <- "Column data"
reddim$ColorByColData <- "condition"

rowdata <- rowDataPlotDefaults(sce, 5)
rowdata$YAxis <- "edgeR.conditiond4Tcf__chir.conditiond4Tcf__unstim.mlog10PValue"
rowdata$XAxis <- "Row data"
rowdata$XAxisRowData <- "edgeR.conditiond4Tcf__chir.conditiond4Tcf__unstim.logFC"

rowstat <- rowStatTableDefaults(sce, 5)
rowstat[['SelectByPlot']] <- c("Row data plot 1", "---", "---", "---", "---")
rowstat[['Selected']] <- c(1L, 2238L, 1L, 1L, 1L)

featassay <- featAssayPlotDefaults(sce, 5)
featassay$XAxis <- "Column data"
featassay$XAxisColData <- "names"
featassay$Assay <- 4L
featassay$PointSize <- 4
featassay[['YAxisRowTable']] <- c("Row statistics table 2", "---", "---", "---", "---")

source("scripts/custom_iSEE_panels.R")
gtf <- prepareGtf("E-MTAB-7029/reference/Homo_sapiens.GRCh38.95.gtf")
saveRDS(gtf, file = "E-MTAB-7029/reference/Homo_sapiens.GRCh38.95.gtf.rds")

cdp <- customDataPlotDefaults(sce, 2)
cdp$Function <- c("customGviz")
cdp$Arguments <- c("bigwig_files E-MTAB-7029/output/STARbigwig/Q10-Chir-1_R1_Aligned.sortedByCoord.out.bw,E-MTAB-7029/output/STARbigwig/Q10-Chir-2_R1_Aligned.sortedByCoord.out.bw,E-MTAB-7029/output/STARbigwig/Q10-Chir-3_R1_Aligned.sortedByCoord.out.bw,E-MTAB-7029/output/STARbigwig/Q10-unstim-1_R1_Aligned.sortedByCoord.out.bw,E-MTAB-7029/output/STARbigwig/Q10-unstim-2_R1_Aligned.sortedByCoord.out.bw,E-MTAB-7029/output/STARbigwig/Q10-unstim-3_R1_Aligned.sortedByCoord.out.bw\nbigwig_names Q10-Chir-1_R1,Q10-Chir-2_R1,Q10-Chir-3_R1,Q10-unstim-1_R1,Q10-unstim-2_R1,Q10-unstim-3_R1\nbigwig_condition d4Tcf__chir,d4Tcf__chir,d4Tcf__chir,d4Tcf__unstim,d4Tcf__unstim,d4Tcf__unstim\ngranges Homo_sapiens.GRCh38.95.gtf.rds\nchr 1\nstart 6.1e6\nend 6.2e6\nshowgene HMOX1")

iSEE(sce, 
     redDimArgs = reddim,
     rowDataArgs = rowdata,
     rowStatArgs = rowstat, 
     featAssayArgs = featassay, 
     customDataArgs = cdp, 
     customDataFun = list(customGviz = customGviz),
     initialPanels = DataFrame(
         Name = c("Reduced dimension plot 1", "Custom data plot 1",
                  "Row data plot 1", "Row statistics table 1", 
                  "Feature assay plot 1", "Row statistics table 2"),
         Width = c(4, 8, 3, 3, 3, 3)
     ))

