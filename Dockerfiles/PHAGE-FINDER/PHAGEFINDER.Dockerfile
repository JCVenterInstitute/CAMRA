FROM --platform=linux/amd64 debian:rc-buggy
LABEL maintainer="Daniella Matute <dmatute@jcvi.org>" 

RUN apt-get update && \
    apt-get install -y wget libssl-dev && \
    apt-get install -y perl cpanminus && \