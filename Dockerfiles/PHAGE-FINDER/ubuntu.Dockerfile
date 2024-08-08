FROM --platform=linux/x86_64 ubuntu:jammy

LABEL author="Daniella Matute" \
    description="A heuristic computer program that identifies prophage regions in completed bacterial genomes." \
    maintainer="dmatute@jcvi.org"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update && \
    apt upgrade -y && \
    apt install -y build-essential wget libssl-dev perl alien cpanminus libxml2 linux-headers-generic blender libgraphics-colorobject-perl libgomp1 python3 python3-pip infernal-doc autoconf git  musl dietlibc-dev && \
    cpanm Math::Round less  Term::ReadKey  Graphics::ColorNames::WWW XML::Simple module

WORKDIR /opt

# PhageFinder
RUN git clone https://github.com/DanyMatute/PhageFinder.git 

# BLAST
ENV BLASTPLUS_VERSION=2.15.0
    # For aarch64-linux architecture
    # RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-aarch64-linux.tar.gz && \
    #     tar -xvf ncbi-blast-2.16.0+-aarch64-linux.tar.gz && \
    #     rm ncbi-blast-2.16.0+-aarch64-linux.tar.gz 
    # ENV PATH="/opt/ncbi-blast-2.16.0+/bin/:{$PATH}"

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-1.x86_64.rpm && \
    alien -k ncbi-blast-2.16.0+-1.x86_64.rpm && \
    dpkg -i ncbi-blast_2.16.0+-1_amd64.deb 

# HMMER 
RUN pip install hmmer==3.4.0.0

# tRNAscan-SE
RUN wget --no-check-certificate http://trna.ucsc.edu/software/trnascan-se-2.0.12.tar.gz  && \ 
    tar -xvf trnascan-se-2.0.12.tar.gz && \
    rm trnascan-se-2.0.12.tar.gz && \
    cd /opt/tRNAscan-SE-2.0 && ./configure --prefix /usr/local && make && make install 
ENV PATH="/usr/local/:$PATH"


# Infernal - instructions from https://github.com/EddyRivasLab/infernal
RUN wget http://eddylab.org/software/infernal/infernal.tar.gz && \     
    tar -zxf infernal.tar.gz && \
    rm infernal.tar.gz && \ 
    cd infernal-1.1.5 && ./configure --prefix /usr/local && make && make check && make install &&\
    cd easel && make install
    
# XGRAPH
RUN wget https://www.xgraph.org/linux/xgraph_4.38_linux64.tar.gz && \ 
    tar -xvf xgraph_4.38_linux64.tar.gz && \
    rm xgraph_4.38_linux64.tar.gz

    
# ARAGRON
RUN mkdir Aragron && cd Aragron && \
    wget https://github.com/TheSEED/aragorn/raw/master/README && \
    wget https://github.com/TheSEED/aragorn/raw/master/aragorn.1 &&\
    wget https://github.com/TheSEED/aragorn/raw/master/aragorn1.2.36.c &&\
    gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.36.c

# libpng-code 
RUN git clone https://git.code.sf.net/p/libpng/code libpng-code && \
    cd libpng-code &&  git checkout origin/libpng15 &&\
    ./configure && make check && make install && ldconfig /usr/local/lib 

WORKDIR /data

