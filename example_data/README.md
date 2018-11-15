## A small RNA-seq example data set

This repository contains a small RNA-seq example data set that may be suitable, e.g., for teaching or testing purposes. The original data files come from the study 

> Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker RM, Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri Jr RA, Tantisira KG, Weiss ST, Lu Q: [RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells.](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099625) PLoS ONE 9(6): e99625. https://doi.org/10.1371/journal.pone.0099625 (2014).

(GEO accession number [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)), in which RNA-Seq was used to characterize the human airway smooth muscle transcriptome at baseline and under asthma treatment conditions. This example data set includes four of the samples from this data set (SRR1039508, SRR1039509, SRR1039512 and SRR1039513), representing Dexamethasone treated and untreated samples from two cell lines (N61311 and N052611). The FASTQ files have been subsetted to include only reads aligning within the first 10M bases of chromosome 1 (of the GRCh38 reference genome). 

In addition to the FASTQ files with the reads, we provide reference annotation files from two sources: Ensembl (release GRCh38.93) and Gencode (v28). For each annotation source, we include a fasta file with the genome sequence, a gtf file with the corresponding gene annotation, and one or more fasta files with transcript sequences. All files are subsetted to include only features from the first 10M bases of chromosome 1. 

### Gencode reference files
The Gencode reference files were downloaded from [https://www.gencodegenes.org/releases/current.html](https://www.gencodegenes.org/releases/current.html). 

- reference/Gencode28/GRCh38.primary_assembly.genome.1.1.10M.fa (genome sequence)
- reference/Gencode28/gencode.v28.transcripts.1.1.10M.fa.gz (transcript sequences)
- reference/Gencode28/gencode.v28.annotation.1.1.10M.gtf (gene annotation)

### Ensembl reference files
The Ensembl reference files were downloaded from [https://www.ensembl.org/info/data/ftp/index.html](https://www.ensembl.org/info/data/ftp/index.html).

- reference/Ensembl.GRCh38.93/Homo_sapiens.GRCh38.dna.chromosome.1.1.10M.fa (genome sequence)
- reference/Ensembl.GRCh38.93/Homo_sapiens.GRCh38.cdna.all.1.1.10M.fa.gz (cDNA transcript sequences)
- reference/Ensembl.GRCh38.93/Homo_sapiens.GRCh38.ncrna.1.1.10M.fa.gz (ncRNA transcript sequences)
- reference/Ensembl.GRCh38.93/Homo_sapiens.GRCh38.93.1.1.10M.gtf (gene annotation)
