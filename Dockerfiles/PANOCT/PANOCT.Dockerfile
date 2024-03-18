FROM --platform=linux/amd64 debian:rc-buggy
LABEL maintainer="Daniella Matute <dmatute@jcvi.org>" dockerfile-date-updated="Jan 11 2024" tool="Panoct" description="Pan-genome ortholog clustering tool (PanOCT) is a tool for pan-genomic analysis of closely related prokaryotic species or strains. PanOCT uses conserved gene neighborhood information to separate recently diverged paralogs into orthologous clusters where homology-only clustering methods cannot."


RUN apt-get update && \
    apt-get install -y perl cpanminus && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN cpanm Data::Dumper
RUN cpanm Cwd
RUN apt-get update && apt-get -y install build-essential wget

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz &&\
    tar -xzf ncbi-blast-2.15.0+-x64-linux.tar.gz &&\
    rm -dr ncbi-blast-2.15.0+-x64-linux.tar.gz

COPY panoct_v3.23.tar.gz /panoct_v3.23.tar.gz
RUN tar -xvzf panoct_v3.23.tar.gz 
RUN rm -f panoct_v3.23.tar.gz 

COPY README.md /README.md
