FROM --platform=linux/amd64 debian:stable-slim


RUN apt update && \
    apt upgrade -y && \
    apt -y install wget


# Install CONDA
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/Miniconda3-latest-Linux-x86_64.sh && \ 
    chmod u+x ~/Miniconda3-latest-Linux-x86_64.sh && \
    ~/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda && \
    . /opt/miniconda/bin/activate
ENV PATH="/opt/miniconda/bin:$PATH"

# Install Platon
RUN conda install -c conda-forge -c bioconda -c defaults platon

# Install DB
RUN wget https://zenodo.org/record/4066768/files/db.tar.gz
RUN tar -xzf db.tar.gz
RUN rm db.tar.gz

WORKDIR /data

