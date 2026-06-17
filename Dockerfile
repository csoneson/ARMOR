FROM continuumio/miniconda3

WORKDIR /app

# Prepare environment
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install -y -c conda-forge mamba=0.23.3
RUN mamba install -y -c bioconda snakemake=6.10.0
