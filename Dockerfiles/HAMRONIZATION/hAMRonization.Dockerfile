FROM --platform=linux/amd64 debian:rc-buggy
LABEL maintainer "Daniella Matute <dmatute@jcvi.org>" 
LABEL base.image="debian:rc-buggy"
LABEL software="hAMRonization"
LABEL software.version="1.1.4"
LABEL dockerfile.updated="Apr 29 2024"
LABEL description="This repo contains the hAMRonization module and CLI parser tools combine the outputs of 18 (as of 2022-09-25) disparate antimicrobial resistance gene detection tools into a single unified format."


RUN apt update 
RUN apt install -y git python3 python3-pip wget 
RUN rm -rf /usr/lib/python3.11/EXTERNALLY-MANAGED
#Install hAMRonization
RUN wget https://github.com/DanyMatute/hAMRonization/archive/refs/tags/v1.1.5.tar.gz && \ 
pip install v1.1.5.tar.gz



WORKDIR /data