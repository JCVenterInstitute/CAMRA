FROM --platform=linux/amd64 debian:rc-buggy
LABEL maintainer="Daniella Matute <dmatute@jcvi.org>" dockerfile-date-updated="Apr 23 2024" tool="hAMRonization" description="This repo contains the hAMRonization module and CLI parser tools combine the outputs of 18 (as of 2022-09-25) disparate antimicrobial resistance gene detection tools into a single unified format."


RUN apt-get update 
RUN apt-get install -y wget libssl-dev
RUN apt install -y git python3 python3-pip python3-pandas

WORKDIR /hAMRonization_workingdir

