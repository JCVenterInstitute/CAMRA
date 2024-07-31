FROM --platform=linux/amd64 debian:stable-slim
LABEL maintainer="Daniella Matute <dmatute@jcvi.org>" 

#Install hmmer2 infernal perl and python muscle
RUN apt update && \
    apt upgrade -y && \
    apt install -y wget libssl-dev perl cpanminus build-essential libxml2 linux-headers-generic blender libgraphics-colorobject-perl libgomp1 python3 python3-pip infernal-doc infernal autoconf hmmer2 git build-essential musl dietlibc-dev && \
    cpanm Math::Round less  Term::ReadKey  Graphics::ColorNames::WWW XML::Simple module


# # PhageFinder
WORKDIR /opt

RUN git clone https://github.com/DanyMatute/PhageFinder.git 

# BLAST
ENV BLASTPLUS_VERSION=2.15.0
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLASTPLUS_VERSION}/ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz && \
    tar -xvf ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz && \
    rm ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz
ENV PATH="/opt/ncbi-blast-2.15.0+/bin/:{$PATH}"

# HMMER 
# RUN wget http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz && \
# tar zxf hmmer-3.3.2.tar.gz && rm hmmer-3.3.2.tar.gz &&\
# cd hmmer-3.3.2&&\
# ./configure --prefix /usr/local && make && make check && make install 


# tRNAscan-SE
# RUN wget --no-check-certificate http://trna.ucsc.edu/software/trnascan-se-2.0.12.tar.gz  && \ 
#     tar -xvf trnascan-se-2.0.12.tar.gz && \
#     rm trnascan-se-2.0.12.tar.gz #&& \
#     cd /opt/tRNAscan-SE-2.0 && ./configure && make && make install &&\ 
#     cd easel && make install
# tRNAscan-SE according to install file
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

WORKDIR /data

