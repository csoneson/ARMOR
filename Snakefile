configfile: "config.yaml"

import pandas as pd
samples = pd.read_table(config["metatxt"])

## Sanitize output directory
import re
def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

outputdir = getpath(config["output"])
FASTQdir = getpath(config["FASTQ"])

print(outputdir)
print(FASTQdir)

## Define the conda environment for all rules using R
if config["useCondaR"] == True:
	Renv = "envs/environment_R.yaml"
else:
	Renv = "envs/environment.yaml"

print(Renv)

Rbin = config["Rbin"]

## ------------------------------------------------------------------------------------ ##
## Target definitions
## ------------------------------------------------------------------------------------ ##
## Run all analyses
## Add outputdir + "outputR/DRIMSeq_dtu.rds" if desired
rule all:
	input:
		outputdir + "MultiQC/multiqc_report.html",
		outputdir + "outputR/edgeR_dge.rds",
		outputdir + "outputR/shiny_results_list.rds",
		outputdir + "outputR/shiny_results_sce.rds",
		outputdir + "outputR/shiny_results_list_edgeR.rds",
		outputdir + "outputR/shiny_results_sce_edgeR.rds"
		
## Install tximeta (or other packages from a repository needed by the user)		
rule pkginstall:
	input:
		script = "scripts/install_pkgs.R"
	output:
	    outputdir + "Rout/pkginstall_state.txt"
	priority: 50
	conda:
		Renv
	log:
		outputdir + "Rout/install_pkgs.Rout"
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args outtxt='{output}' " {input.script} {log}'''

## FastQC on original (untrimmed) files
rule runfastqc:
	input:
		expand(outputdir + "FastQC/{sample}_R1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_R2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

## Trimming and FastQC on trimmed files
rule runtrimming:
	input:
		expand(outputdir + "FastQC/{sample}_R1_val_1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_R2_val_2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_trimmed_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist())

## Salmon quantification
rule runsalmonquant:
	input:
		expand(outputdir + "salmon/{sample}/quant.sf", sample = samples.names.values.tolist())

## STAR alignment
rule runstar:
	input:
		expand(outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.names.values.tolist())

## List all the packages that were used by the R analyses
rule listpackages:
	log:
		outputdir + "Rout/list_packages.Rout"
	params:
		Routdir = outputdir + "Rout",
		outtxt = outputdir + "R_package_versions.txt",
		script = "scripts/list_packages.R"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args Routdir='{params.Routdir}' outtxt='{params.outtxt}'" {params.script} {log}'''

## Print the versions of all software packages
rule softwareversions:
	conda:
		"envs/environment.yaml"
	shell:
		"salmon --version; trim_galore --version; cutadapt --version; "
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
		outputdir + "logs/salmon_index.log"
	params:
		salmonk = config["salmonk"],
		salmonoutdir = config["salmonindex"],
		anno = config["annotation"]
	conda:
		"envs/environment.yaml"
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
		script = "scripts/generate_linkedTxome.R",
		install = outputdir + "Rout/pkginstall_state.txt"
	log:
		outputdir + "Rout/generate_linkedTxome.Rout"
	output:
		config["salmonindex"] + ".json"
	params:
		flag = config["annotation"],
		organism = config["organism"],
		release = str(config["release"]),
		build = config["build"]
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args transcriptfasta='{input.txome}' salmonidx='{input.salmonidx}' gtf='{input.gtf}' annotation='{params.flag}' organism='{params.organism}' release='{params.release}' build='{params.build}' output='{output}'" {input.script} {log}'''

## Generate STAR index
rule starindex:
	input:
		genome = config["genome"],
		gtf = config["gtf"]
	output:
		config["STARindex"] + "/SA",
		config["STARindex"] + "/chrNameLength.txt"
	log:
		outputdir + "logs/STAR_index.log"
	params:
		STARindex = config["STARindex"],
		readlength = config["readlength"]
	conda:
		"envs/environment.yaml"
	threads: 
		config["ncores"]
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
		fastq = FASTQdir + "{sample}.fastq.gz"
	output:
		outputdir + "FastQC/{sample}_fastqc.zip"
	params:
		FastQC = outputdir + "FastQC"
	log:
		outputdir + "logs/fastqc_{sample}.log"
	conda:
		"envs/environment.yaml"
	threads: 
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## FastQC, trimmed reads
rule fastqc2:
	input:
		fastq = outputdir + "FASTQtrimmed/{sample}.fq.gz"
	output:
		outputdir + "FastQC/{sample}_fastqc.zip"
	params:
		FastQC = outputdir + "FastQC"
	log:
		outputdir + "logs/fastqc_trimmed_{sample}.log"
	conda:
		"envs/environment.yaml"
	threads: 
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## MultiQC
rule multiqc:
	input:
		expand(outputdir + "FastQC/{sample}_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_R1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_R2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_trimmed_fastqc.zip", sample = samples.names[samples.type == 'SE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_R1_val_1_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FastQC/{sample}_R2_val_2_fastqc.zip", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz", sample = samples.names[samples.type == 'SE'].values.tolist()),
		expand(outputdir + "FASTQtrimmed/{sample}_R1_val_1.fq.gz", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "FASTQtrimmed/{sample}_R2_val_2.fq.gz", sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(outputdir + "salmon/{sample}/quant.sf", sample = samples.names.values.tolist()),
		expand(outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.names.values.tolist())
	output:
		outputdir + "MultiQC/multiqc_report.html"
	params:
		MultiQCdir = outputdir + "MultiQC",
		FastQCdir = outputdir + "FastQC",
		FASTQtrimmeddir = outputdir + "FASTQtrimmed",
		salmondir = outputdir + "salmon",
		STARdir = outputdir + "STAR" 
	log:
		outputdir + "logs/multiqc.log"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.FastQCdir} {params.FASTQtrimmeddir} {params.salmondir} {params.STARdir} -f -o {params.MultiQCdir}"


## ------------------------------------------------------------------------------------ ##
## Adapter trimming
## ------------------------------------------------------------------------------------ ##
# TrimGalore!
rule trimgaloreSE:
	input:
		fastq = FASTQdir + "{sample}.fastq.gz"
	output:
		outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz"
	params:
		FASTQtrimmeddir = outputdir + "FASTQtrimmed"
	log:
		outputdir + "logs/trimgalore_{sample}.log"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq}"

rule trimgalorePE:
	input:
		fastq1 = FASTQdir + "{sample}_R1.fastq.gz",
		fastq2 = FASTQdir + "{sample}_R2.fastq.gz"
	output:
		outputdir + "FASTQtrimmed/{sample}_R1_val_1.fq.gz",
		outputdir + "FASTQtrimmed/{sample}_R2_val_2.fq.gz"
	params:
		FASTQtrimmeddir = outputdir + "FASTQtrimmed"
	log:
		outputdir + "logs/trimgalore_{sample}.log"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt "
		"--paired {input.fastq1} {input.fastq2}"

## ------------------------------------------------------------------------------------ ##
## Salmon abundance estimation
## ------------------------------------------------------------------------------------ ##
# Estimate abundances with Salmon
rule salmonSE:
	input:
		index = config["salmonindex"] + "/hash.bin",
		fastq = outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz"
	output:
		outputdir + "salmon/{sample}/quant.sf"
	log:
		outputdir + "logs/salmon_{sample}.log"
	threads: 
		config["ncores"]
	params:
		salmonindex = config["salmonindex"],
		fldMean = config["fldMean"],
		fldSD = config["fldSD"],
		salmondir = outputdir + "salmon"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'Salmon version:\n' > {log}; salmon --version >> {log}; "
		"salmon quant -i {params.salmonindex} -l A -r {input.fastq} "
		"-o {params.salmondir}/{wildcards.sample} --seqBias --gcBias "
		"--fldMean {params.fldMean} --fldSD {params.fldSD} -p {threads}"

rule salmonPE:
	input:
		index = config["salmonindex"] + "/hash.bin",
		fastq1 = outputdir + "FASTQtrimmed/{sample}_R1_val_1.fq.gz",
		fastq2 = outputdir + "FASTQtrimmed/{sample}_R2_val_2.fq.gz"
	output:
		outputdir + "salmon/{sample}/quant.sf"
	log:
		outputdir + "logs/salmon_{sample}.log"
	threads: 
		config["ncores"]
	params:
		salmonindex = config["salmonindex"],
		fldMean = config["fldMean"],
		fldSD = config["fldSD"],
		salmondir = outputdir + "salmon"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'Salmon version:\n' > {log}; salmon --version >> {log}; "
		"salmon quant -i {params.salmonindex} -l A -1 {input.fastq1} -2 {input.fastq2} "
		"-o {params.salmondir}/{wildcards.sample} --seqBias --gcBias "
		"--fldMean {params.fldMean} --fldSD {params.fldSD} -p {threads}"

## ------------------------------------------------------------------------------------ ##
## STAR mapping
## ------------------------------------------------------------------------------------ ##
## Genome mapping with STAR
rule starSE:
	input:
		index = config["STARindex"] + "/SA",
		fastq = outputdir + "FASTQtrimmed/{sample}_trimmed.fq.gz"
	output:
		outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	threads: 
		config["ncores"]
	log:
		outputdir + "logs/STAR_{sample}.log"
	params:
		STARindex = config["STARindex"],
		STARdir = outputdir + "STAR"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c"

rule starPE:
	input:
		index = config["STARindex"] + "/SA",
		fastq1 = outputdir + "FASTQtrimmed/{sample}_R1_val_1.fq.gz",
		fastq2 = outputdir + "FASTQtrimmed/{sample}_R2_val_2.fq.gz"
	output:
		outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	threads: 
		config["ncores"]
	log:
		outputdir + "logs/STAR_{sample}.log"
	params:
		STARindex = config["STARindex"],
		STARdir = outputdir + "STAR"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq1} {input.fastq2} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c"

## Index bam files
rule staridx:
	input:
		bam = outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
	log:
		outputdir + "logs/samtools_index_{sample}.log"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"

## Convert BAM files to bigWig
rule bigwig:
	input:
		bam = outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
		chrl = config["STARindex"] + "/chrNameLength.txt"
	output:
		outputdir + "STARbigwig/{sample}_Aligned.sortedByCoord.out.bw"
	params:
		STARbigwigdir = outputdir + "STARbigwig"
	log:
		outputdir + "logs/bigwig_{sample}.log"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -ibam {input.bam} -bg | sort -k1,1 -k2,2n > "
		"{params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph; "
		"bedGraphToBigWig {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph "
		"{input.chrl} {output}; rm -f {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph"

## ------------------------------------------------------------------------------------ ##
## Transcript quantification
## ------------------------------------------------------------------------------------ ##
## tximeta
rule tximeta:
	input:
	    outputdir + "Rout/pkginstall_state.txt",
		expand(outputdir + "salmon/{sample}/quant.sf", sample = samples.names.values.tolist()),
		metatxt = config["metatxt"],
		salmonidx = config["salmonindex"] + "/hash.bin",
		json = config["salmonindex"] + ".json",
		script = "scripts/run_tximeta.R"
	output:
		outputdir + "outputR/tximeta_se.rds"
	log:
		outputdir + "Rout/tximeta_se.Rout"
	params:
		salmondir = outputdir + "salmon"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args salmondir='{params.salmondir}' json='{input.json}' metafile='{input.metatxt}' outrds='{output}'" {input.script} {log}'''


## ------------------------------------------------------------------------------------ ##
## Differential expression
## ------------------------------------------------------------------------------------ ##
## edgeR
rule edgeR:
	input:
	    outputdir + "Rout/pkginstall_state.txt",
		rds = outputdir + "outputR/tximeta_se.rds",
		script = "scripts/run_dge_edgeR.R"
	output:
		outputdir + "outputR/edgeR_dge.rds"
	log:
		outputdir + "Rout/run_dge_edgeR.Rout"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' outrds='{output}'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## Differential transcript usage
## ------------------------------------------------------------------------------------ ##
## DRIMSeq
rule DRIMSeq:
	input:
	    outputdir + "Rout/pkginstall_state.txt",
		rds = outputdir + "outputR/tximeta_se.rds",
		script = "scripts/run_dtu_drimseq.R"
	output:
		outputdir + "outputR/DRIMSeq_dtu.rds"
	log:
		outputdir + "Rout/run_dtu_drimseq.Rout"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' outrds='{output}'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## Shiny app
## ------------------------------------------------------------------------------------ ##
rule shiny:
	input:
	    outputdir + "Rout/pkginstall_state.txt",
		expand(outputdir + "STARbigwig/{sample}_Aligned.sortedByCoord.out.bw", sample = samples.names.values.tolist()),
		rds = outputdir + "outputR/edgeR_dge.rds",
		metatxt = config["metatxt"],
		gtf = config["gtf"],
		script = "scripts/prepare_results_for_shiny.R"
	log: 
		outputdir + "Rout/shiny_results.Rout"
	output:
		outList = outputdir + "outputR/shiny_results_list.rds",
		outSCE = outputdir + "outputR/shiny_results_sce.rds"
	params:
		groupvar = config["groupvar"],
		bigwigdir = "STARbigwig"
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args edgerres='{input.rds}' groupvar='{params.groupvar}' gtffile='{input.gtf}' metafile='{input.metatxt}' bigwigdir='{params.bigwigdir}' outList='{output.outList}' outSCE='{output.outSCE}'" {input.script} {log}'''


rule shinyedgeR:
	input:
	    outputdir + "Rout/pkginstall_state.txt",
		rds = outputdir + "outputR/edgeR_dge.rds",
		metatxt = config["metatxt"],
		script = "scripts/prepare_results_for_shiny.R"
	log:
		outputdir + "Rout/shiny_results_edgeR.Rout"
	output:
		outList = outputdir + "outputR/shiny_results_list_edgeR.rds",
		outSCE = outputdir + "outputR/shiny_results_sce_edgeR.rds"
	params:
		groupvar = config["groupvar"]
	conda:
		Renv
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args edgerres='{input.rds}' groupvar='{params.groupvar}' gtffile=NULL metafile='{input.metatxt}' bigwigdir=NULL outList='{output.outList}' outSCE='{output.outSCE}'" {input.script} {log}'''
