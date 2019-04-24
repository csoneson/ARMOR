# run this script within a reasonably-named current working directory
# (e.g., "E-MTAB-7029" as is done here)

# create directories
dir.create("reference")
dir.create("FASTQ")

# download files and make metadata file
md <- read.table("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7029/E-MTAB-7029.sdrf.txt", 
                 header=TRUE, sep="\t")[,c(1,28,31)]
dummy <- apply(md, 1, function(u) download.file(u[["Comment.FASTQ_URI."]], 
                                                file.path("FASTQ",u[["Scan.Name"]])))

metadata <- data.frame(names = gsub(".fastq.gz","",md$Scan.Name), 
                       type="SE", 
                       condition=gsub("[1-3]$","",md$Source.Name))
write.table(metadata, "metadata.txt", sep="\t", quote=FALSE, row.names=FALSE)

## https://www.ebi.ac.uk/ena/browse/read-download#downloading_files_aspera
## alternative to download.file()
# stub <- "~/.aspera/connect/bin/ascp -QT -P33001"
# stub <- paste0(stub, " -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh")
# stub <- paste0(stub, " era-fasp@fasp.sra.ebi.ac.uk:")
# dummy <- apply(md, 1, function(u) {
#  cmd <- paste0(stub, gsub("ftp://ftp.sra.ebi.ac.uk","",u[["Comment.FASTQ_URI."]]))
#  cmd <- paste0(cmd, " ", file.path("FASTQ",u[["Scan.Name"]]))
#  print(cmd)
#  system(cmd)
#})


# download reference files
ref_files <- c("ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
               "ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
               "ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz")
dummy <- lapply(ref_files, function(u) download.file(u, file.path("reference",basename(u))))

# make sure the genome FA and GTF are unzipped
for(i in c(1,3))
  system(paste0("gunzip reference/", basename(ref_files[i])))
