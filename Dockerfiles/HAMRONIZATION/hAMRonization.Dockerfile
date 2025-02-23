FROM debian:rc-buggy
LABEL maintainer="Daniella Matute <dmatute@jcvi.org>" 
LABEL base.image="debian:rc-buggy"
LABEL software="hAMRonization"
LABEL software.version="1.1.4"
LABEL dockerfile.updated="Apr 29 2024"
LABEL description="This repo contains the hAMRonization module and CLI parser tools combine the outputs of 18 (as of 2022-09-25) disparate antimicrobial resistance gene detection tools into a single unified format."


RUN apt update 
RUN apt install -y git python3 python3-pip wget python3-biopython python3-tabulate python3-setuptools
RUN rm -rf /usr/lib/python3.11/EXTERNALLY-MANAGED

RUN cd /opt && \
    wget https://github.com/dhimmel/obonet/archive/refs/tags/v1.1.0.tar.gz && \
    tar -xvzf v1.1.0.tar.gz && \
    rm v1.1.0.tar.gz && \
    cd obonet-1.1.0 && \
    python3 setup.py install

# Install hAMRonization
RUN cd /opt && \ 
    mkdir AMR_Term_Consolidation && \
    wget https://github.com/DanyMatute/hAMRonization/archive/refs/tags/v1.2.0.tar.gz && \ 
    tar -xzf v1.2.0.tar.gz  && \ 
    rm v1.2.0.tar.gz && \ 
    cd hAMRonization-1.2.0 && \ 
    apt install -y python3-setuptools && \ 
    python3 setup.py install 


# AMR Term Consolidation 
COPY AMR_Term_Consolidation /opt/AMR_Term_Consolidation

WORKDIR /data