FROM continuumio/miniconda3
MAINTAINER Katharina Hembach (katharina.hembach@uzh.ch)

## SHELL ["/bin/bash"]

## system requirements
RUN apt-get update \
    && apt-get install -y htop

## snakemake installation
RUN conda update -n base conda
RUN conda install -c bioconda -c conda-forge snakemake

## copy the environment file to the image and create the conda environment
COPY envs/environment.yaml /tmp/
RUN conda env create -n rnaseqworkflow --file tmp/environment.yaml 

##  start the conda environment in each container
ENV PATH /opt/conda/envs/rnaseqworkflow/bin:$PATH
## alternative: RUN echo "source activate rnaseqworkflow" > ~/.bashrc

## clean up
RUN apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /opt/conda/pkgs/*

CMD [ "/bin/bash" ] 
