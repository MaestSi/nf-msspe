FROM continuumio/miniconda3

########### set variables
ENV DEBIAN_FRONTEND noninteractive

######### dependencies
RUN apt-get update -qq \
    && apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    nano \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

############################################################ install nf-msspe
WORKDIR /home/

RUN conda config --add channels r && \
conda config --add channels anaconda && \
conda config --add channels conda-forge && \
conda config --add channels bioconda

RUN conda install mafft bioconductor-biostrings
