ARG RESFINDER_VER="4.5.0"
FROM ubuntu:jammy 
LABEL maintainer="Daniella Matute <dmatute@jcvi.org>" dockerfile-date-updated="Jan 11 2024" tool="GGRaSP" description="GGRaSP (Gaussian Genome Representative Selector with Prioritization), a R-package and associated executable Rscript program that generates a list of prioritized representative genomes from either supervised or unsupervised clustering of related genomes."


RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    procps \
    git \
    ncbi-blast+ \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-dev \
    gcc \
    make \
    libz-dev \
    unzip && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

#RUN pipx install resfinder==4.5.0 --include-deps && mkdir /data
RUN wget https://bitbucket.org/genomicepidemiology/resfinder/get/842f6193de0d6921c41761ea2caeb0e62d5b8037.tar.gz && \
tar -xvzf 842f6193de0d6921c41761ea2caeb0e62d5b8037.tar.gz && \
rm 842f6193de0d6921c41761ea2caeb0e62d5b8037.tar.gz && \
mv genomicepidemiology-resfinder-842f6193de0d/ resfinder4.5.0 && \ 
apt update && apt install -y python3-biopython python3-cgecore python3-tabulate 
#python3-cgelib


RUN git clone --branch 1.4.14 --depth 1 https://bitbucket.org/genomicepidemiology/kma.git && \
cd kma && \
make && \
mv -v kma* /usr/local/bin/

RUN mkdir /pointfinder_db && \
git clone --branch 4.1.0 --depth 1 https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git /pointfinder_db && \
cd /pointfinder_db && \
python3 INSTALL.py /usr/local/bin/kma_index non_interactive

# resfinder database install
RUN mkdir /resfinder_db && \
git clone --branch 2.3.1 --depth 1 https://git@bitbucket.org/genomicepidemiology/resfinder_db.git /resfinder_db && \
cd /resfinder_db && \
python3 INSTALL.py /usr/local/bin/kma_index non_interactive

# disinfinder database install
RUN mkdir /disinfinder_db && \
git clone --branch 2.0.1 --depth 1 https://git@bitbucket.org/genomicepidemiology/disinfinder_db.git /disinfinder_db && \
cd /disinfinder_db && \
python3 INSTALL.py /usr/local/bin/kma_index non_interactive

ENV PATH="${PATH}" \
# specifying database locations
CGE_RESFINDER_RESGENE_DB="/resfinder_db" \
CGE_RESFINDER_RESPOINT_DB="/pointfinder_db" \
CGE_DISINFINDER_DB="/disinfinder_db"


# setting a janky alias for everyone that uses the "latest" tag 
# RUN echo -e '#!/bin/bash\npython3 -m resfinder "$@"' > /usr/bin/run_resfinder.py && \
# chmod +x /usr/bin/run_resfinder.py