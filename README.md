# cool-name RNA-seq workflow

cool-name is a [Snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html), aimed at performing a typical RNA-seq workflow in a reproducible, automated, and partially contained manner. It is implemented such that alternative or similar analysis can be added or removed.  

cool-name consists of a `Snakefile`, a [`conda`](https://conda.io/docs/) environment file (`envs/environment.yaml`) a configuration file (`config.yaml`) and a set of `R` scripts, to perform quality control, preprocessing and differential expression analysis of RNA-seq data. The output can be combined with the [`iSEE`](https://github.com/csoneson/iSEE) `R` package to generate a `shiny` application for browsing and sharing the results.  

This initiative was proposed by [Charlotte Soneson](https://github.com/csoneson).
Current contributors include:

- [Ruizhu Huang](https://github.com/fionarhuang)
- [Katharina Hembach](https://github.com/khembach)
- [Stephany Orjuela](https://github.com/sorjuela)
- [Mark D. Robinson](https://github.com/markrobinsonuzh)
- [Charlotte Soneson](https://github.com/csoneson)

## Workflow graph
![DAG](dag_nice3.png)  
Blue circles are rules run in `R`, orange circles from software called as shell commands. Dashed lines and light-colored circles are optional rules, controlled in `config.yaml`


## Using the RNA-seq workflow
To use the RNA-seq workflow on your own data, follow the steps outlined in the [wiki](https://github.com/csoneson/iSEE/wiki).
