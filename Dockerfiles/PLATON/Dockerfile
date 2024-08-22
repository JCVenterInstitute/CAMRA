FROM --platform=linux/amd64 debian:stable-slim


RUN apt update && \
    apt upgrade -y && \
    apt -y install wget ncbi-blast+ build-essential


# Install DB
WORKDIR /var/lib

RUN mkdir platon_db && cd platon_db && \
    wget https://zenodo.org/record/4066768/files/db.tar.gz && \
    tar -xzf db.tar.gz && \
    rm db.tar.gz
WORKDIR /

# Install CONDA
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/Miniconda3-latest-Linux-x86_64.sh && \ 
    chmod u+x ~/Miniconda3-latest-Linux-x86_64.sh && \
    ~/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda && \
    . /opt/miniconda/bin/activate
ENV PATH="/opt/miniconda/bin:$PATH"

WORKDIR /opt


#5
RUN wget https://files.pythonhosted.org/packages/b9/70/7d09c2b4144c5bb328e3b4a74ee2f3ba63bf95e625a60bd9cf5b839d37bf/cb-platon-1.7.tar.gz && \ 
tar -xvf cb-platon-1.7.tar.gz && rm cb-platon-1.7.tar.gz &&  cd cb-platon-1.7 && python setup.py install

#Prodigal (2.6.3)
RUN conda install -y -c conda-forge -c bioconda prodigal=2.6.3 

#Diamon (2.0.6)
RUN conda install -y -c conda-forge -c bioconda diamond=2.0.14 

#Blast+ (2.10.1) 
# version 2.12.0+ds-3 is the only version available in apt and it works with Platon. Installed above. 

#Mummer
RUN wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz && \
    tar -xzf mummer-4.0.0beta2.tar.gz && \    
    rm mummer-4.0.0beta2.tar.gz && \
    cd mummer-4.0.0beta2 && \
    ./configure && \
    make  && \
    make install && \
    ldconfig

#HMMER (3.3.1)
RUN conda install --solver=classic  -y -c conda-forge -c bioconda hmmer=3.3.1

#INFERNAL (1.1.4)
RUN conda install --solver=classic -y -c conda-forge -c bioconda infernal=1.1.4

WORKDIR /data

