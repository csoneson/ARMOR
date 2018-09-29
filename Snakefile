configfile: "config.yaml"

import pandas as pd
samples = pd.read_table(config["metatxt"])

## ------------------------------------------------------------------------------------ ##
## Target definitions
## ------------------------------------------------------------------------------------ ##
## Run all analyses
## Add "output/DRIMSeq_dtu.rds" if desired
rule all:
	input:
		config["output"]+"/MultiQC/multiqc_report.html",
		config["output"]+"/outputR/edgeR_dge.rds",
		config["output"]+"/outputR/shiny_results_list.rds",
		config["output"]+"/outputR/shiny_results_sce.rds",
		config["output"]+"/outputR/shiny_results_list_edgeR.rds",
		config["output"]+"/outputR/shiny_results_sce_edgeR.rds"

## FastQC on original (untrimmed) files
rule runfastqc:
	input:
 		expand(config["output"]+"/FastQC/{sample}_R1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
 		expand(config["output"]+"/FastQC/{sample}_R2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(config["output"]+"/FastQC/{sample}_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

## Trimming and FastQC on trimmed files
rule runtrimming:
	input:
 		expand(config["output"]+"/FastQC/{sample}_R1_val_1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
 		expand(config["output"]+"/FastQC/{sample}_R2_val_2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(config["output"]+"/FastQC/{sample}_trimmed_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

## Salmon quantification
rule runsalmonquant:
	input:
		expand(config["output"]+"/salmon/{sample}/quant.sf", sample = samples.names.values.tolist())

## STAR alignment
rule runstar:
	input:
		expand(config["output"]+"/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.names.values.tolist())

## List all the packages that were used by the R analyses
rule listpackages:
	log:
		config["output"]+"/Rout/list_packages.Rout"
	params:
		Routdir = "Rout",
		outtxt = "R_package_versions.txt",
		script = "scripts/list_packages.R"
	shell:
		'''R CMD BATCH --no-restore --no-save "--args Routdir='{params.Routdir}' outtxt='{params.outtxt}'" {input.script} {log}'''

## Print the versions of all software packages
rule softwareversions:
	shell:
		"R --version; salmon --version; trim_galore --version; cutadapt --version; "
		"fastqc --version; STAR --version; samtools --version; multiqc --version; "
		"bedtools --version"

## ------------------------------------------------------------------------------------ ##
## Reference preparation
## ------------------------------------------------------------------------------------ ##
## Generate Salmon index from merged cDNA and ncRNA files
rule salmonindex:
	input:
		txome = config["txome"]
	output:
		config["salmonindex"] + "/hash.bin"
	log:
		config["output"]+"/logs/salmon_index.log"
	params:
		salmonk = config["salmonk"],
		salmonoutdir = config["salmonindex"],
		anno =  config["annotation"]
	shell:
	  """
	  if [ {params.anno} == "Gencode" ]; then
      echo 'Salmon version:\n' > {log}; salmon --version >> {log};
  	  salmon index -t {input.txome} -k {params.salmonk} -i {params.salmonoutdir} --gencode --type quasi

    else
  	  echo 'Salmon version:\n' > {log}; salmon --version >> {log};
      salmon index -t {input.txome} -k {params.salmonk} -i {params.salmonoutdir} --type quasi
    fi
    """

## Generate linkedTxome mapping
rule linkedTxome:
	input:
		txome = config["txome"],
		gtf = config["gtf"],
		salmonidx = config["salmonindex"] + "/hash.bin",
		script = "scripts/generate_linkedTxome.R"
	log:
		config["output"]+"/Rout/generate_linkedTxome.Rout"
	output:
	  config["salmonindex"] + ".json"
	params:
		flag = config["annotation"],
		organism = config["organism"],
		release = str(config["release"]),
		build = config["build"]
	shell:
		'''R CMD BATCH --no-restore --no-save "--args transcriptfasta='{input.txome}' salmonidx='{input.salmonidx}' gtf='{input.gtf}' annotation='{params.flag}' organism='{params.organism}' release='{params.release}' build='{params.build}' output='{output}'" {input.script} {log}'''

## Generate STAR index
rule starindex:
	input:
		genome = config["genome"],
		gtf = config["gtf"]
	output:
		config["STARindex"] + "/SA",
		config["STARindex"] + "/chrNameLength.txt"
	log:
		config["output"]+"/logs/STAR_index.log"
	params:
		STARindex = config["STARindex"],
		readlength = config["readlength"]
	threads: config["ncores"]
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.STARindex} "
		"--genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang {params.readlength}"

## ------------------------------------------------------------------------------------ ##
## Quality control
## ------------------------------------------------------------------------------------ ##
## FastQC, original reads
rule fastqc:
	input:
		fastq = config["FASTQ"]+"/{sample}.fastq.gz"
	output:
		config["output"]+"/FastQC/{sample}_fastqc.zip"
	params:
	    FastQC = config["output"]+"/FastQC"
	log:
		config["output"]+"/logs/fastqc_{sample}.log"
	threads: config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## FastQC, trimmed reads
rule fastqc2:
	input:
		fastq = config["output"] +"/FASTQtrimmed/{sample}.fq.gz"
	output:
		config["output"]+"/FastQC/{sample}_fastqc.zip"
	params:
	    FastQC = config["output"]+"/FastQC"
	log:
		config["output"]+"/logs/fastqc_trimmed_{sample}.log"
	threads: config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## MultiQC
rule multiqc:
	input:
		expand(config["output"]+"/FastQC/{sample}_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist()),
		expand(config["output"]+"/FastQC/{sample}_R1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(config["output"]+"/FastQC/{sample}_R2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(config["output"]+"/FastQC/{sample}_trimmed_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist()),
		expand(config["output"]+"/FastQC/{sample}_R1_val_1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(config["output"]+"/FastQC/{sample}_R2_val_2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(config["output"] +"/FASTQtrimmed/{sample}_trimmed.fq.gz", sample = samples.names[samples.type == 'SE'].values.tolist()),
		expand(config["output"] +"/FASTQtrimmed/{sample}_R1_val_1.fq.gz", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(config["output"] +"/FASTQtrimmed/{sample}_R2_val_2.fq.gz", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(config["output"]+"/salmon/{sample}/quant.sf", sample = samples.names.values.tolist()),
		expand(config["output"]+"/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.names.values.tolist())
	output:
		config["output"]+"/MultiQC/multiqc_report.html"
	params:
	    MultiQC = config["output"]+"/MultiQC",
	    FastQC = config["output"]+"/FastQC",
	    FASTQtrimmed = config["output"]+"/FASTQtrimmed",
	    salmon = config["output"]+"/salmon",
	    STAR = config["output"]+"/STAR" 
	log:
		config["output"]+"/logs/multiqc.log"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.FastQC} {params.FASTQtrimmed} {params.salmon} {params.STAR} -f -o {params.MultiQC}"


## ------------------------------------------------------------------------------------ ##
## Adapter trimming
## ------------------------------------------------------------------------------------ ##
# TrimGalore!
rule trimgaloreSE:
	input:
		fastq = config["FASTQ"]+"/{sample}.fastq.gz"
	output:
		config["output"] +"/FASTQtrimmed/{sample}_trimmed.fq.gz"
	params:
	    FASTQtrimmed = config["output"] + "/FASTQtrimmed"
	log:
		config["output"]+"/logs/trimgalore_{sample}.log"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmed} --path_to_cutadapt cutadapt {input.fastq}"

rule trimgalorePE:
	input:
		fastq1 = config["FASTQ"]+"/{sample}_R1.fastq.gz",
		fastq2 = config["FASTQ"]+"/{sample}_R2.fastq.gz"
	output:
	    config["output"] +"/FASTQtrimmed/{sample}_R1_val_1.fq.gz",
		config["output"] +"/FASTQtrimmed/{sample}_R2_val_2.fq.gz"
	params:
	    FASTQtrimmed = config["output"] + "/FASTQtrimmed"
	log:
		config["output"]+"/logs/trimgalore_{sample}.log"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmed} --path_to_cutadapt cutadapt "
		"--paired {input.fastq1} {input.fastq2}"

## ------------------------------------------------------------------------------------ ##
## Salmon abundance estimation
## ------------------------------------------------------------------------------------ ##
# Estimate abundances with Salmon
rule salmonSE:
	input:
		index = config["salmonindex"] + "/hash.bin",
		fastq = config["output"] +"/FASTQtrimmed/{sample}_trimmed.fq.gz"
	output:
		config["output"]+"/salmon/{sample}/quant.sf"
	log:
		config["output"]+"/logs/salmon_{sample}.log"
	threads: config["ncores"]
	params:
		salmonindex = config["salmonindex"],
		fldMean = config["fldMean"],
		fldSD = config["fldSD"],
		salmon = config["output"]+"/salmon"
	shell:
		"echo 'Salmon version:\n' > {log}; salmon --version >> {log}; "
		"salmon quant -i {params.salmonindex} -l A -r {input.fastq} "
		"-o {params.salmon}/{wildcards.sample} --seqBias --gcBias "
		"--fldMean {params.fldMean} --fldSD {params.fldSD} -p {threads}"

rule salmonPE:
	input:
		index = config["salmonindex"] + "/hash.bin",
		fastq1 = config["output"] +"/FASTQtrimmed/{sample}_R1_val_1.fq.gz",
		fastq2 = config["output"] +"/FASTQtrimmed/{sample}_R2_val_2.fq.gz"
	output:
		config["output"]+"/salmon/{sample}/quant.sf"
	log:
		config["output"]+"/logs/salmon_{sample}.log"
	threads: config["ncores"]
	params:
		salmonindex = config["salmonindex"],
		fldMean = config["fldMean"],
		fldSD = config["fldSD"],
		salmon = config["output"]+"/salmon"
	shell:
		"echo 'Salmon version:\n' > {log}; salmon --version >> {log}; "
		"salmon quant -i {params.salmonindex} -l A -1 {input.fastq1} -2 {input.fastq2} "
		"-o {params.salmon}/{wildcards.sample} --seqBias --gcBias "
		"--fldMean {params.fldMean} --fldSD {params.fldSD} -p {threads}"

## ------------------------------------------------------------------------------------ ##
## STAR mapping
## ------------------------------------------------------------------------------------ ##
## Genome mapping with STAR
rule starSE:
	input:
		index = config["STARindex"] + "/SA",
		fastq = config["output"] +"/FASTQtrimmed/{sample}_trimmed.fq.gz"
	output:
		config["output"]+"/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	threads: config["ncores"]
	log:
		config["output"]+"/logs/STAR_{sample}.log"
	params:
		STARindex = config["STARindex"],
		STAR = config["output"] +"/STAR"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq} "
		"--runThreadN {threads} --outFileNamePrefix {params.STAR}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c"

rule starPE:
	input:
		index = config["STARindex"] + "/SA",
		fastq1 = config["output"] +"/FASTQtrimmed/{sample}_R1_val_1.fq.gz",
		fastq2 = config["output"] +"/FASTQtrimmed/{sample}_R2_val_2.fq.gz"
	output:
		config["output"]+"/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	threads: config["ncores"]
	log:
		config["output"]+"/logs/STAR_{sample}.log"
	params:
		STARindex = config["STARindex"],
		STAR = config["output"] +"/STAR"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq1} {input.fastq2} "
		"--runThreadN {threads} --outFileNamePrefix {params.STAR}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c"

## Index bam files
rule staridx:
	input:
		bam = config["output"]+"/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		config["output"]+"/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
	log:
		config["output"]+"/logs/samtools_index_{sample}.log"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"

## Convert BAM files to bigWig
rule bigwig:
	input:
		bam = config["output"]+"/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
		chrl = config["STARindex"] + "/chrNameLength.txt"
	output:
		config["output"]+"/STARbigwig/{sample}_Aligned.sortedByCoord.out.bw"
	params:
	    STARbigwig = config["output"]+"/STARbigwig"
	log:
		config["output"]+"/logs/bigwig_{sample}.log"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -ibam {input.bam} -bg | sort -k1,1 -k2,2n > "
		"{params.STARbigwig}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph; "
		"bedGraphToBigWig {params.STARbigwig}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph "
		"{input.chrl} {output}; rm -f {params.STARbigwig}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph"

## ------------------------------------------------------------------------------------ ##
## Differential expression
## ------------------------------------------------------------------------------------ ##
## edgeR
rule edgeR:
	input:
		expand(config["output"]+"/salmon/{sample}/quant.sf", sample = samples.names.values.tolist()),
		metatxt = config["metatxt"],
		salmonidx = config["salmonindex"] + "/hash.bin",
		json = config["salmonindex"] + ".json",
		script = "scripts/run_dge_edgeR.R"
	output:
		config["output"]+"/outputR/edgeR_dge.rds"
	log:
		config["output"]+"/Rout/run_dge_edgeR.Rout"
	params:
		salmon = config["output"]+"/salmon",
	shell:
		'''R CMD BATCH --no-restore --no-save "--args salmondir='{params.salmon}' json='{input.json}' metafile='{input.metatxt}' outrds='{output}'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## Differential transcript usage
## ------------------------------------------------------------------------------------ ##
## DRIMSeq
rule DRIMSeq:
	input:
		expand(config["output"]+"/salmon/{sample}/quant.sf", sample = samples.names.values.tolist()),
		metatxt = config["metatxt"],
		script = "scripts/run_dtu_drimseq.R"
	output:
		config["output"]+"/outputR/DRIMSeq_dtu.rds"
	log:
		config["output"]+"/Rout/run_dtu_drimseq.Rout"
	params:
		salmon = config["output"]+"/salmon",
	shell:
		'''R CMD BATCH --no-restore --no-save "--args salmondir='{params.salmon}' metafile='{input.metatxt}' outrds='{output}'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## Shiny app
## ------------------------------------------------------------------------------------ ##
rule shiny:
	input:
		expand(config["output"]+"/STARbigwig/{sample}_Aligned.sortedByCoord.out.bw", sample = samples.names.values.tolist()),
		rds = config["output"]+"/outputR/edgeR_dge.rds",
		metatxt = config["metatxt"],
		gtf = config["gtf"],
		script = "scripts/prepare_results_for_shiny.R"
	log: config["output"]+"/Rout/shiny_results.Rout"
	output:
		outList = config["output"]+"/outputR/shiny_results_list.rds",
		outSCE = config["output"]+"/outputR/shiny_results_sce.rds"
	params:
		groupvar = config["groupvar"],
		bigwigdir = "STARbigwig"
	shell:
		'''R CMD BATCH --no-restore --no-save "--args edgerres='{input.rds}' groupvar='{params.groupvar}' gtffile='{input.gtf}' metafile='{input.metatxt}' bigwigdir='{params.bigwigdir}' outList='{output.outList}' outSCE='{output.outSCE}'" {input.script} {log}'''


rule shinyedgeR:
	input:
		rds = config["output"]+"/outputR/edgeR_dge.rds",
		metatxt = config["metatxt"],
		script = "scripts/prepare_results_for_shiny.R"
	log:
		config["output"]+"/Rout/shiny_results_edgeR.Rout"
	output:
		outList = config["output"]+"/outputR/shiny_results_list_edgeR.rds",
		outSCE = config["output"]+"/outputR/shiny_results_sce_edgeR.rds"
	params:
		groupvar = config["groupvar"]
	shell:
		'''R CMD BATCH --no-restore --no-save "--args edgerres='{input.rds}' groupvar='{params.groupvar}' gtffile=NULL metafile='{input.metatxt}' bigwigdir=NULL outList='{output.outList}' outSCE='{output.outSCE}'" {input.script} {log}'''
