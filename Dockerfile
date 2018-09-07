FROM continuumio/miniconda3
MAINTAINER Katharina Hembach (katharina.hembach@uzh.ch)

## SHELL ["/bin/bash"]

## system requirements
RUN apt-get update \
    && apt-get install -y htop


## snakemake installation
RUN conda update -n base conda
RUN conda install -c bioconda -c conda-forge snakemake


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





##### alternative: copy the whole github repo
## clone the rnaseqworkflow from github
## RUN git clone https://github.com/csoneson/rnaseqworkflow.git

## make the conda environment with all required software
## TODO: or do we want to have an image with all software installed without a conda environment?
##WORKDIR rnaseqworkflow
##RUN conda env create -n rnaseqworkflow --file envs/environment.yaml 








### Ideas:

## TODO: do not clone the whole repo but just copy the environment file to the image?
## the user has to start the container from the github repo and mount it to be able to run the pipeline?
## COPY envs/environment.yaml /tmp/
## RUN conda env create -n rnaseqworkflow --file tmp/environment.yaml 


## TODO: how do we get the metadata.txt, the user defined config.yaml, the FASTQ files and the reference annotation into the docker container?


## mount the directory and then call docker with
## docker run --name rnaseqworkflow -it --mount type=bind,source="$(pwd)",target=/rnaseqworkflow khembach/rnaseqworkflow bin/bash
## start snakemake in the interactive session --> preferred way!!



## run container in detached mode
## docker run --name rnaseqworkflow -dit --mount type=bind,source="$(pwd)",target=/rnaseqworkflow khembach/rnaseqworkflow
## execute snakemake in running container
## docker exec rnaseqworkflow snakemake runfastqc
## --> you need to make sure that you are in the directory where you mounted the data!!



## commit to new branch
git checkout master
git pull
git checkout <your_branch>
git merge master

