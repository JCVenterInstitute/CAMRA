FROM --platform=linux/amd64 debian:stable-slim
LABEL maintainer="Daniella Matute <dmatute@jcvi.org>" 

RUN apt update && \
    apt upgrade -y && \
    apt install -y wget libssl-dev perl cpanminus libxml2 libgomp1 python3 python3-pip infernal-doc autoconf git build-essential musl dietlibc-dev && \
    cpanm Math::Round less 

ENV BIN=/bin
ENV OPT=/opt

# # PhageFinder
WORKDIR /opt

####4####
#COPY ../../../PhageFinder /opt/phage_finder
####3####
RUN git clone https://github.com/DanyMatute/PhageFinder.git
####2####
#COPY ./phage_finder_v2.5_4docker /opt/phage_finder_v2.5_4
####1####
# COPY ./phage_finder_v2.5_4docker.tar.gz /opt

# RUN tar -xvf phage_finder_v2.5_4docker.tar.gz && \
#     rm phage_finder_v2.5_4docker.tar.gz

#ENV PATH="/opt/phage_finder_v2.5_4docker/bin:{$PATH}"ß


# BLAST
ENV BLASTPLUS_VERSION=2.15.0
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLASTPLUS_VERSION}/ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz && \
    tar -xvf ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz && \
    rm ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz
ENV PATH="/opt/ncbi-blast-2.15.0+/bin/:{$PATH}"


# HMMER
RUN apt-get update && (apt-get install -t buster-backports -y hmmer || apt-get install -y hmmer) && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

# tRNAscan-SE
RUN wget --no-check-certificate http://trna.ucsc.edu/software/trnascan-se-2.0.12.tar.gz  && \ 
    tar -xvf trnascan-se-2.0.12.tar.gz && \
    rm trnascan-se-2.0.12.tar.gz #&& \
    cd /opt/tRNAscan-SE-2.0 && ./configure && make && make install


# XGRAPH
RUN wget https://www.xgraph.org/linux/xgraph_4.38_linux64.tar.gz && \ 
    tar -xvf xgraph_4.38_linux64.tar.gz && \
    rm xgraph_4.38_linux64.tar.gz

    
# ARAGRON
RUN mkdir Aragron && cd Aragron && \
    wget https://github.com/TheSEED/aragorn/raw/master/README && \
    wget https://github.com/TheSEED/aragorn/raw/master/aragorn.1 &&\
    wget https://github.com/TheSEED/aragorn/raw/master/aragorn1.2.36.c

# SEQSTAT
RUN git clone --recursive https://github.com/clwgg/seqstats &&\
    cd seqstats && make

#EASEL (for esl-seqstat)    
RUN git clone https://github.com/EddyRivasLab/easel && \ 
    cd easel && autoconf && ./configure && make && make check \
    wget /opt/easel 

RUN wget https://uclibc.org/downloads/uClibc-0.9.33.2.tar.xz &&\
    tar -xvf uClibc-0.9.33.2.tar.xz &&\
    rm uClibc-0.9.33.2.tar.xz 


WORKDIR /data

