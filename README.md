# ARMOR RNA-seq workflow

**ARMOR** (Automated Reproducible MOdular Rna) is a [Snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html), aimed at performing a typical RNA-seq workflow in a reproducible, automated, and partially contained manner. It is implemented such that alternative or similar analysis can be added or removed. 

ARMOR consists of a `Snakefile`, a [`conda`](https://conda.io/docs/) environment file (`envs/environment.yaml`) a configuration file (`config.yaml`) and a set of `R` scripts, to perform quality control, preprocessing and differential expression analysis of RNA-seq data. The output can be combined with the [`iSEE`](https://github.com/csoneson/iSEE) `R` package to generate a `shiny` application for browsing and sharing the results.  

By default, the pipeline performs the steps shown in the [diagram](dag_nice3.png) below. However, if there is a specific step you do not want to run (e.g `STAR` alignment), it can be easily "turned off" in the `config.yaml`. 

*Advanced use*: If you prefer other software to run one of the outlined steps (e.g. `DEXSeq` over `edgeR`, or `kallisto` over `Salmon`), you can use the software of your preference provided you have your own script(s), and change some lines within the `Snakefile`. If you think your "custom rule" might be of use to a broader audience, let us know by opening an issue.


## Using the RNA-seq workflow
To use the RNA-seq workflow on your own data, follow the steps outlined in the [wiki](https://github.com/csoneson/rnaseqworkflow/wiki).

## Workflow graph
![DAG](dag_nice3.png)  
Blue circles are rules run in `R`, orange circles from software called as shell commands. Dashed lines and light-colored circles are optional rules, controlled in `config.yaml`

## Contributors
This initiative was proposed by [Charlotte Soneson](https://github.com/csoneson).
Current contributors include:

- [Ruizhu Huang](https://github.com/fionarhuang)
- [Katharina Hembach](https://github.com/khembach)
- [Stephany Orjuela](https://github.com/sorjuela)
- [Mark D. Robinson](https://github.com/markrobinsonuzh)
- [Charlotte Soneson](https://github.com/csoneson)