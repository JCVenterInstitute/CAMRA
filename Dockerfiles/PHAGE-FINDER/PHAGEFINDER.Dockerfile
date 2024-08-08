FROM --platform=linux/amd64 debian:stable-slim
LABEL maintainer="Daniella Matute <dmatute@jcvi.org>" 

RUN apt update && \
    apt upgrade -y && \
    apt install -y wget libssl-dev perl cpanminus libxml2 && \
    cpanm Math::Round less

ENV BIN /bin
ENV OPT /opt

# # PhageFinder
WORKDIR /opt

COPY ./phage_finder_v2.5_4docker.tar.gz /opt
RUN tar -xvf phage_finder_v2.5_4docker.tar.gz && \
    rm phage_finder_v2.5_4docker.tar.gz

# BLAST
ENV BLASTPLUS_VERSION 2.15.0
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLASTPLUS_VERSION}/ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz && \
    tar -xvf ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz && \
    rm ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz

# # # AB-BLAST
# COPY ./AB-BLAST /bin/AB-BLAST
# WORKDIR /bin/AB-BLAST
# RUN mkdir ~/.config && mkdir ~/.config/ab-blast && mv license.xml ~/.config/ab-blast && \
#     tar -xvf ab-blast-20200317-linux-x64.tar.gz  && \
#     rm ab-blast-20200317-linux-x64.tar.gz && \
#     chmod 600 ~/.config/ab-blast/license.xml

# HMMER
RUN apt-get update && (apt-get install -t buster-backports -y hmmer || apt-get install -y hmmer) && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

# tRNAscan-SE
RUN wget --no-check-certificate http://trna.ucsc.edu/software/trnascan-se-2.0.12.tar.gz  && \ 
    tar -xvf trnascan-se-2.0.12.tar.gz && \
    rm trnascan-se-2.0.12.tar.gz


# XGRAPH
RUN wget https://www.xgraph.org/linux/xgraph_4.38_linux64.tar.gz && \ 
    tar -xvf xgraph_4.38_linux64.tar.gz && \
    rm xgraph_4.38_linux64.tar.gz

# ARAGRON
RUN mkdir Aragron && cd Aragron && \
    wget https://github.com/TheSEED/aragorn/raw/master/README && \
    wget https://github.com/TheSEED/aragorn/raw/master/aragorn.1 &&\
    wget https://github.com/TheSEED/aragorn/raw/master/aragorn1.2.36.c


WORKDIR /data

