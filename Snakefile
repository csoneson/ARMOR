## Configuration file
#from pathlib import Path
#from os import path
import os
if len(config) == 0:
  if os.path.isfile("./config.yaml"):
    configfile: "./config.yaml"
  else:
    sys.exit("Make sure there is a config.yaml file in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")

## Make sure that all expected variables from the config file are in the config dictionary
configvars = ['annotation', 'organism', 'build', 'release', 'txome', 'genome', 'gtf', 'salmonindex', 'salmonk', 'STARindex', 'readlength', 'fldMean', 'fldSD', 'metatxt', 'design', 'contrast', 'genesets', 'ncores', 'FASTQ', 'fqext1', 'fqext2', 'fqsuffix', 'output', 'useCondaR', 'Rbin', 'run_trimming', 'run_STAR', 'run_DRIMSeq', 'run_camera']
for k in configvars:
	if k not in config:
		config[k] = None

## If any of the file paths is missing, replace it with ""
def sanitizefile(str):
	if str is None:
		str = ''
	return str

config['txome'] = sanitizefile(config['txome'])
config['gtf'] = sanitizefile(config['gtf'])
config['genome'] = sanitizefile(config['genome'])
config['STARindex'] = sanitizefile(config['STARindex'])
config['salmonindex'] = sanitizefile(config['salmonindex'])
config['metatxt'] = sanitizefile(config['metatxt'])

## Read metadata
if not os.path.isfile(config["metatxt"]):
  sys.exit("Metadata file " + config["metatxt"] + " does not exist.")

import pandas as pd
samples = pd.read_csv(config["metatxt"], sep='\t')

if not set(['names','type']).issubset(samples.columns):
  sys.exit("Make sure 'names' and 'type' are columns in " + config["metatxt"])


## Sanitize provided input and output directories
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

## Define the conda environment for all rules using R
if config["useCondaR"] == True:
	Renv = "envs/environment_R.yaml"
else:
	Renv = "envs/environment.yaml"

## Define the R binary
Rbin = config["Rbin"]

## ------------------------------------------------------------------------------------ ##
## Target definitions
## ------------------------------------------------------------------------------------ ##
## Run all analyses
rule all:
	input:
		os.path.join(outputdir, "MultiQC", "multiqc_report.html"),
		os.path.join(outputdir, "outputR", "shiny_sce.rds")

rule setup:
	input:
		f"{outputdir}Rout/pkginstall_state.txt",
		f"{outputdir}Rout/softwareversions.done"

## Install R packages
rule pkginstall:
	input:
		script = "scripts/install_pkgs.R"
	output:
	  	f"{outputdir}Rout/pkginstall_state.txt"
	params:
		flag = config["annotation"],
		ncores = config["ncores"],
		organism = config["organism"],
		Rbin = Rbin
	priority:
		50
	conda:
		Renv
	log:
		f"{outputdir}Rout/install_pkgs.Rout"
	benchmark:
	  	f"{outputdir}benchmarks/install_pkgs.txt"
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args outtxt='{output}' ncores='{params.ncores}' annotation='{params.flag}' organism='{params.organism}'" {input.script} {log}'''

## FastQC on original (untrimmed) files
rule runfastqc:
	input:
		expand(os.path.join(outputdir,"FastQC/{sample}_", str(config["fqext1"]), "_fastqc.zip"), sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(os.path.join(outputdir,"FastQC/{sample}_", str(config["fqext2"]), "_fastqc.zip"), sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(os.path.join(outputdir,"FastQC/{sample}_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist())

## Trimming and FastQC on trimmed files
rule runtrimming:
	input:
		expand(os.path.join(outputdir, "FastQC/{sample}_", str(config["fqext1"]), "_val_1_fastqc.zip"), sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(os.path.join(outputdir, "FastQC/{sample}_", str(config["fqext2"]), "_val_2_fastqc.zip"), sample = samples.names[samples.type == 'PE'].values.tolist()),
		expand(os.path.join(outputdir, "FastQC","{sample}_trimmed_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist())

## Salmon quantification
rule runsalmonquant:
	input:
		expand(os.path.join(outputdir, "salmon/{sample}/quant.sf"), sample = samples.names.values.tolist())

## STAR alignment
rule runstar:
	input:
		expand(os.path.join(outputdir, "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"), sample = samples.names.values.tolist()),
		expand(os.path.join(outputdir, "STARbigwig/{sample}_Aligned.sortedByCoord.out.bw"), sample = samples.names.values.tolist())

## List all the packages that were used by the R analyses
rule listpackages:
	log:
		f"{outputdir}Rout/list_packages.Rout"
	params:
		Routdir = f"{outputdir}Rout",
		outtxt = f"{outputdir}R_package_versions.txt",
		script = "scripts/list_packages.R",
		Rbin = Rbin
	conda:
		Renv
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args Routdir='{params.Routdir}' outtxt='{params.outtxt}'" {params.script} {log}'''

## Print the versions of all software packages
rule softwareversions:
	output:
		touch(f"{outputdir}Rout/softwareversions.done")
	log:
		f"{outputdir}logs/softversions.log"
	conda:
		"envs/environment.yaml"
	shell:
		"echo -n 'ARMOR version ' && cat version; "
		"salmon --version; trim_galore --version; "
		"echo -n 'cutadapt ' && cutadapt --version; "
		"fastqc --version; STAR --version; samtools --version; multiqc --version; "
		"bedtools --version"

## ------------------------------------------------------------------------------------ ##
## Reference preparation
## ------------------------------------------------------------------------------------ ##
## Generate Salmon index from merged cDNA and ncRNA files
rule salmonindex:
	input:
		txome = config["txome"],
		salmonoutdir = config["salmonindex"]
	output:
		os.path.join(config["salmonindex"], "versionInfo.json")
	log:
		f"{outputdir}logs/salmon_index.log"
	benchmark:
		f"{outputdir}benchmarks/salmon_index.txt"
	params:
		#salmonk = config["salmonk"],
		anno = config["annotation"],
		salmonextraparams = config["additional_salmon_index"]
	conda:
		"envs/environment.yaml"
	shell:
	  """
	  if [ {params.anno} == "Gencode" ]; then
      echo 'Salmon version:\n' > {log}; salmon --version >> {log};
  	  salmon index -t {input.txome} -i {params.salmonoutdir} --gencode {params.salmonextraparams}

    else
  	  echo 'Salmon version:\n' > {log}; salmon --version >> {log};
      salmon index -t {input.txome} -i {input.salmonoutdir} {params.salmonextraparams}
    fi
    """

## Generate linkedtxome mapping
rule linkedtxome:
	input:
		txome = config["txome"],
		gtf = config["gtf"],
		salmonidx = os.path.join(config["salmonindex"], "versionInfo.json"),
		script = "scripts/generate_linkedtxome.R",
		install = f"{outputdir}Rout/pkginstall_state.txt"
	log:
		f"{outputdir}Rout/generate_linkedtxome.Rout"
	benchmark:
		f"{outputdir}benchmarks/generate_linkedtxome.txt"
	output:
		config["salmonindex"] + ".json"
	params:
		flag = config["annotation"],
		organism = config["organism"],
		release = str(config["release"]),
		build = config["build"],
		Rbin = Rbin
	conda:
		Renv
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args transcriptfasta='{input.txome}' salmonidx='{input.salmonidx}' gtf='{input.gtf}' annotation='{params.flag}' organism='{params.organism}' release='{params.release}' build='{params.build}' output='{output}'" {input.script} {log}'''

## Generate STAR index
rule starindex:
	input:
		genome = config["genome"],
		gtf = config["gtf"],
		STARindex = config["STARindex"]
	output:
		os.path.join(config["STARindex"], "SA"),
		os.path.join(config["STARindex"], "chrNameLength.txt")
	log:
		f"{outputdir}logs/STAR_index.log"
	benchmark:
		f"{outputdir}benchmarks/STAR_index.txt"
	params:
		readlength = config["readlength"],
		starextraparams = config["additional_star_index"]
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {input.STARindex} "
		"--genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang {params.readlength} "
		"{params.starextraparams}"

## ------------------------------------------------------------------------------------ ##
## Quality control
## ------------------------------------------------------------------------------------ ##
## FastQC, original reads
rule fastqc:
	input:
		fastq = os.path.join(FASTQdir, "{sample}." + str(config["fqsuffix"]) + ".gz")
	output:
		os.path.join(outputdir, "FastQC/{sample}_fastqc.zip")
	params:
		FastQC = f"{outputdir}FastQC"
	log:
		os.path.join(outputdir, "logs/fastqc_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks/fastqc_{sample}.txt")
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"

## FastQC, trimmed reads
rule fastqctrimmed:
	input:
		fastq = os.path.join(outputdir, "FASTQtrimmed/{sample}.fq.gz")
	output:
		os.path.join(outputdir, "FastQC/{sample}_fastqc.zip")
	params:
		FastQC = os.path.join(outputdir,"FastQC")
	log:
		os.path.join(outputdir, "logs/fastqc_trimmed_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks/fastqc_trimmed_{sample}.txt")
	conda:
		"envs/environment.yaml"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o {params.FastQC} -t {threads} {input.fastq}"



# The config.yaml files determines which steps should be performed
def multiqc_input(wildcards):
	input = []
	input.extend(expand(os.path.join(outputdir, "FastQC/{sample}_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist()))
	input.extend(expand(os.path.join(outputdir, "FastQC/{sample}_" + str(config["fqext1"]) + "_fastqc.zip"), sample = samples.names[samples.type == 'PE'].values.tolist()))
	input.extend(expand(os.path.join(outputdir, "FastQC/{sample}_" + str(config["fqext2"]) + "_fastqc.zip"), sample = samples.names[samples.type == 'PE'].values.tolist()))
	input.extend(expand(os.path.join(outputdir, "salmon/{sample}/quant.sf"), sample = samples.names.values.tolist()))
	if config["run_trimming"]:
		input.extend(expand(os.path.join(outputdir, "FASTQtrimmed/{sample}_trimmed.fq.gz"), sample = samples.names[samples.type == 'SE'].values.tolist()))
		input.extend(expand(os.path.join(outputdir, "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz"), sample = samples.names[samples.type == 'PE'].values.tolist()))
		input.extend(expand(os.path.join(outputdir, "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz"), sample = samples.names[samples.type == 'PE'].values.tolist()))
		input.extend(expand(os.path.join(outputdir, "FastQC/{sample}_trimmed_fastqc.zip"), sample = samples.names[samples.type == 'SE'].values.tolist()))
		input.extend(expand(os.path.join(outputdir, "FastQC/{sample}_" + str(config["fqext1"]) + "_val_1_fastqc.zip"), sample = samples.names[samples.type == 'PE'].values.tolist()))
		input.extend(expand(os.path.join(outputdir, "FastQC/{sample}_" + str(config["fqext2"]) + "_val_2_fastqc.zip"), sample = samples.names[samples.type == 'PE'].values.tolist()))
	if config["run_STAR"]:
		input.extend(expand(os.path.join(outputdir, "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"), sample = samples.names.values.tolist()))
	return input

## Determine the input directories for MultiQC depending on the config file
def multiqc_params(wildcards):
	param = [outputdir + "FastQC",
	outputdir + "salmon"]
	if config["run_trimming"]:
		param.append(outputdir + "FASTQtrimmed")
	if config["run_STAR"]:
		param.append(outputdir + "STAR")
	return param

## MultiQC
rule multiqc:
	input:
		multiqc_input
	output:
		f"{outputdir}MultiQC/multiqc_report.html"
	params:
		inputdirs = multiqc_params,
		MultiQCdir = f"{outputdir}MultiQC"
	log:
		f"{outputdir}logs/multiqc.log"
	benchmark:
		f"{outputdir}benchmarks/multiqc.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc {params.inputdirs} -f -o {params.MultiQCdir}"


## ------------------------------------------------------------------------------------ ##
## Adapter trimming
## ------------------------------------------------------------------------------------ ##
# TrimGalore!
rule trimgaloreSE:
	input:
		fastq = os.path.join(FASTQdir, "{sample}." + str(config["fqsuffix"]) + ".gz")
	output:
		os.path.join(outputdir, "FASTQtrimmed/{sample}_trimmed.fq.gz")
	params:
		FASTQtrimmeddir = f"{outputdir}FASTQtrimmed"
	log:
		os.path.join(outputdir, "logs/trimgalore_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks/trimgalore_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o {params.FASTQtrimmeddir} --path_to_cutadapt cutadapt {input.fastq}"

rule trimgalorePE:
	input:
		fastq1 = os.path.join(FASTQdir, "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]) + ".gz"),
		fastq2 = os.path.join(FASTQdir, "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]) + ".gz")
	output:
		os.path.join(outputdir, "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz"),
		os.path.join(outputdir, "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz")
	params:
		FASTQtrimmeddir = f"{outputdir}FASTQtrimmed"
	log:
		os.path.join(outputdir, "logs/trimgalore_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks/trimgalore_{sample}.txt")
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
		index = os.path.join(config["salmonindex"], "versionInfo.json"),
		fastq = os.path.join(outputdir, "FASTQtrimmed/{sample}_trimmed.fq.gz") if config["run_trimming"] else os.path.join(FASTQdir, "{sample}." + str(config["fqsuffix"]) + ".gz")
	output:
		outputdir + "salmon/{sample}/quant.sf"
	log:
		os.path.join(outputdir, "logs/salmon_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks/salmon_{sample}.txt")
	threads:
		config["ncores"]
	params:
		salmonindex = config["salmonindex"],
		salmondir = f"{outputdir}salmon",
		salmonextraparams = config["additional_salmon_quant"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'Salmon version:\n' > {log}; salmon --version >> {log}; "
		"salmon quant -i {params.salmonindex} -l A -r {input.fastq} "
		"-o {params.salmondir}/{wildcards.sample} -p {threads} {params.salmonextraparams}"

rule salmonPE:
	input:
		index = os.path.join(config["salmonindex"], "versionInfo.json"),
		fastq1 = os.path.join(outputdir, "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz") if config["run_trimming"] else os.path.join(FASTQdir, "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]) + ".gz"),
		fastq2 = os.path.join(outputdir, "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz") if config["run_trimming"] else os.path.join(FASTQdir, "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]) + ".gz")
	output:
		outputdir + "salmon/{sample}/quant.sf"
	log:
		os.path.join(outputdir, "logs/salmon_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks/salmon_{sample}.txt")
	threads:
		config["ncores"]
	params:
		salmonindex = config["salmonindex"],
		salmondir = f"{outputdir}salmon",
		salmonextraparams = config["additional_salmon_quant"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'Salmon version:\n' > {log}; salmon --version >> {log}; "
		"salmon quant -i {params.salmonindex} -l A -1 {input.fastq1} -2 {input.fastq2} "
		"-o {params.salmondir}/{wildcards.sample} -p {threads} {params.salmonextraparams}"

## ------------------------------------------------------------------------------------ ##
## STAR mapping
## ------------------------------------------------------------------------------------ ##
## Genome mapping with STAR
rule starSE:
	input:
		index = os.path.join(config["STARindex"],"SA"),
		fastq = os.path.join(outputdir, "FASTQtrimmed/{sample}_trimmed.fq.gz") if config["run_trimming"] else os.path.join(FASTQdir, "{sample}." + str(config["fqsuffix"]) + ".gz")
	output:
		outputdir + "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	threads:
		config["ncores"]
	log:
		os.path.join(outputdir, "logs/STAR_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks/STAR_{sample}.txt")
	params:
		STARindex = config["STARindex"],
		STARdir = f"{outputdir}STAR",
		starextraparams = config["additional_star_align"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
		"{params.starextraparams}"

rule starPE:
	input:
		index = os.path.join(config["STARindex"], "SA"),
		fastq1 = os.path.join(outputdir, "FASTQtrimmed/{sample}_" + str(config["fqext1"]) + "_val_1.fq.gz") if config["run_trimming"] else os.path.join(FASTQdir, "{sample}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]) + ".gz"),
		fastq2 = os.path.join(outputdir, "FASTQtrimmed/{sample}_" + str(config["fqext2"]) + "_val_2.fq.gz") if config["run_trimming"] else os.path.join(FASTQdir, "{sample}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]) + ".gz")
	output:
		os.path.join(outputdir, "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam")
	threads:
		config["ncores"]
	log:
		os.path.join(outputdir, "logs/STAR_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks/STAR_{sample}.txt")
	params:
		STARindex = config["STARindex"],
		STARdir = f"{outputdir}STAR",
		starextraparams = config["additional_star_align"]
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq1} {input.fastq2} "
		"--runThreadN {threads} --outFileNamePrefix {params.STARdir}/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
		"{params.starextraparams}"

## Index bam files
rule bamindex:
	input:
		bam = os.path.join(outputdir, "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam")
	output:
		os.path.join(outputdir, "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai")
	log:
		os.path.join(outputdir, "logs/samtools_index_{sample}.log")
	benchmark:
		os.path.join(outputdir, "benchmarks/samtools_index_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"

## Convert BAM files to bigWig
rule bigwig:
	input:
		bam = os.path.join(outputdir, "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"),
		chrl = os.path.join(config["STARindex"], "chrNameLength.txt")
	output:
		os.path.join(outputdir, "STARbigwig/{sample}_Aligned.sortedByCoord.out.bw")
	params:
		STARbigwigdir = f"{outputdir}STARbigwig"
	log:
		outputdir + "logs/bigwig_{sample}.log"
	benchmark:
		os.path.join(outputdir, "benchmarks/bigwig_{sample}.txt")
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'bedtools version:\n' > {log}; bedtools --version >> {log}; "
		"bedtools genomecov -split -ibam {input.bam} -bg | LC_COLLATE=C sort -k1,1 -k2,2n > "
		"{params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph; "
		"bedGraphToBigWig {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph "
		"{input.chrl} {output}; rm -f {params.STARbigwigdir}/{wildcards.sample}_Aligned.sortedByCoord.out.bedGraph"

## ------------------------------------------------------------------------------------ ##
## Transcript quantification
## ------------------------------------------------------------------------------------ ##
## tximeta
rule tximeta:
	input:
	    f"{outputdir}Rout/pkginstall_state.txt",
		expand(os.path.join(outputdir, "salmon/{sample}/quant.sf"), sample = samples.names.values.tolist()),
		metatxt = config["metatxt"],
		salmonidx = os.path.join(config["salmonindex"], "versionInfo.json"),
		json = config["salmonindex"] + ".json",
		script = "scripts/run_tximeta.R"
	output:
		f"{outputdir}outputR/tximeta_se.rds"
	log:
		f"{outputdir}Rout/tximeta_se.Rout"
	benchmark:
		f"{outputdir}benchmarks/tximeta_se.txt"
	params:
		salmondir = f"{outputdir}salmon",
		flag = config["annotation"],
		organism = config["organism"],
		Rbin = Rbin
	conda:
		Renv
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args salmondir='{params.salmondir}' json='{input.json}' metafile='{input.metatxt}' outrds='{output}' annotation='{params.flag}' organism='{params.organism}'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## Input variable check
## ------------------------------------------------------------------------------------ ##
def geneset_param(wildcards):
	if config["run_camera"]:
                gs = config["genesets"].replace(" ", "") if config["genesets"] is not None else "NOTDEFINED"
		return "genesets='" + gs + "'"
	else:
		return ""


## check design matrix and contrasts
rule checkinputs:
    input:
        "config.yaml",
        script = "scripts/check_input.R"
    output:
        f"{outputdir}Rout/check_input.txt"
    log:
        f"{outputdir}Rout/check_input.Rout"
    benchmark:
    	f"{outputdir}benchmarks/check_input.txt"
    params:
        gtf = config["gtf"],
        genome = config["genome"],
        txome = config["txome"],
        fastqdir = config["FASTQ"],
        metatxt = config["metatxt"],
        design = config["design"].replace(" ", "") if config["design"] is not None else "NOTDEFINED",
        contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else "NOTDEFINED",
        annotation = config["annotation"].replace(" ", "") if config["annotation"] is not None else "NOTDEFINED",
        genesets = geneset_param,
        fqsuffix = str(config["fqsuffix"]),
        fqext1 = str(config["fqext1"]),
        fqext2 = str(config["fqext2"]),
        run_camera = str(config["run_camera"]),
        organism = config["organism"],
	Rbin = Rbin
    conda:
	    Renv
    shell:
        '''{params.Rbin} CMD BATCH --no-restore --no-save "--args metafile='{params.metatxt}' design='{params.design}' contrast='{params.contrast}' outFile='{output}' gtf='{params.gtf}' genome='{params.genome}' fastqdir='{params.fastqdir}' fqsuffix='{params.fqsuffix}' fqext1='{params.fqext1}' fqext2='{params.fqext2}' txome='{params.txome}' run_camera='{params.run_camera}' organism='{params.organism}' {params.genesets} annotation='{params.annotation}'" {input.script} {log};
        cat {output}
        '''
       

## ------------------------------------------------------------------------------------ ##
## Differential expression
## ------------------------------------------------------------------------------------ ##
rule edgeR:
	input:
		f"{outputdir}Rout/pkginstall_state.txt",
		rds = f"{outputdir}outputR/tximeta_se.rds",
		script = "scripts/run_render.R",
		template = "scripts/edgeR_dge.Rmd"
	output:
		html = f"{outputdir}outputR/edgeR_dge.html",
		rds = f"{outputdir}outputR/edgeR_dge.rds"
	params:
		directory = f"{outputdir}outputR",
		organism = config["organism"],        
                design = config["design"].replace(" ", "") if config["design"] is not None else "",
                contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else "",
		genesets = geneset_param,
		Rbin = Rbin
	log:
		f"{outputdir}Rout/run_dge_edgeR.Rout"
	benchmark:
		f"{outputdir}benchmarks/run_dge_edgeR.txt"
	conda:
		Renv
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' organism='{params.organism}' design='{params.design}' contrast='{params.contrast}' {params.genesets} rmdtemplate='{input.template}' outputdir='{params.directory}' outputfile='edgeR_dge.html'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## Differential transcript usage
## ------------------------------------------------------------------------------------ ##
## DRIMSeq
rule DRIMSeq:
	input:
	    f"{outputdir}Rout/pkginstall_state.txt",
		rds = f"{outputdir}outputR/edgeR_dge.rds",
		script = "scripts/run_render.R",
		template = "scripts/DRIMSeq_dtu.Rmd"
	output:
		html = f"{outputdir}outputR/DRIMSeq_dtu.html",
		rds = f"{outputdir}outputR/DRIMSeq_dtu.rds"
	params:
		directory = f"{outputdir}outputR",
		organism = config["organism"],
		ncores = config["ncores"],
                design = config["design"].replace(" ", "") if config["design"] is not None else "",
                contrast = config["contrast"].replace(" ", "") if config["contrast"] is not None else "",
		Rbin = Rbin
	log:
		f"{outputdir}Rout/run_dtu_drimseq.Rout"
	benchmark:
		f"{outputdir}benchmarks/run_dtu_drimseq.txt"
	conda:
		Renv
	threads:
		config["ncores"]
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' design='{params.design}' contrast='{params.contrast}' ncores='{params.ncores}' rmdtemplate='{input.template}' outputdir='{params.directory}' outputfile='DRIMSeq_dtu.html'" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## shiny app
## ------------------------------------------------------------------------------------ ##
def shiny_input(wildcards):
	input = [f"{outputdir}Rout/pkginstall_state.txt"]
	if config["run_STAR"]:
		input.extend(expand(os.path.join(outputdir, "STARbigwig/{sample}_Aligned.sortedByCoord.out.bw"), sample = samples.names.values.tolist()))
	return input

def shiny_params(wildcards):
	param = ["outputdir='" + outputdir + "outputR'"]
	if config["run_STAR"]:
		param.append("bigwigdir='" + outputdir + "STARbigwig'")
	return param

## shiny
rule shiny:
	input:
		shiny_input,
		rds = f"{outputdir}outputR/DRIMSeq_dtu.rds" if config["run_DRIMSeq"]
			else f"{outputdir}outputR/edgeR_dge.rds",
		script = "scripts/run_render.R",
		gtf = config["gtf"],
		template = "scripts/prepare_shiny.Rmd"
	output:
		html = f"{outputdir}outputR/prepare_shiny.html",
		rds = f"{outputdir}outputR/shiny_sce.rds"
	params:
		p = shiny_params,
		Rbin = Rbin
	log:
		f"{outputdir}Rout/prepare_shiny.Rout"
	benchmark:
		f"{outputdir}benchmarks/prepare_shiny.txt"
	conda:
		Renv
	shell:
		'''{params.Rbin} CMD BATCH --no-restore --no-save "--args se='{input.rds}' gtffile='{input.gtf}' rmdtemplate='{input.template}' outputfile='prepare_shiny.html' {params.p}" {input.script} {log}'''

## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##
onsuccess:
	print("Success! The Snakemake workflow is completed.")

onerror:
	print("Error! The Snakemake workflow aborted.")
