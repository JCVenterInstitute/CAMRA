ARG RESFINDER_VER="4.5.0"

FROM --platform=linux/amd64 debian:rc-buggy

LABEL maintainer="Daniella Matute <dmatute@jcvi.org>" dockerfile-date-updated="Jan 11 2024" tool="Resfinder" description="ResFinder identifies (BLAST & KMA) acquired antimicrobial resistance genes in total or partial sequenced isolates of bacteria."

ARG RESFINDER_VER
ARG KMA_VER="1.4.14"
ARG RESFINDER_DB_COMMIT_HASH="2.3.1"
ARG POINTFINDER_DB_COMMIT_HASH="4.1.0"
ARG DISINFINDER_DB_COMMIT_HASH="2.0.1"


# Install deprendencies 
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


# Install resfinder 
RUN pip3 install --break-system-packages resfinder==${RESFINDER_VER} && mkdir /data

# Install KMA
RUN git clone --branch ${KMA_VER} --depth 1 https://bitbucket.org/genomicepidemiology/kma.git && \
cd kma && \
make && \
mv -v kma* /usr/local/bin/

# pointfinder db
RUN mkdir /pointfinder_db && \
git clone --branch ${POINTFINDER_DB_COMMIT_HASH} --depth 1 https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git /pointfinder_db && \
cd /pointfinder_db && \
python3 INSTALL.py /usr/local/bin/kma_index non_interactive

# resfinder database install
RUN mkdir /resfinder_db && \
git clone --branch ${RESFINDER_DB_COMMIT_HASH}  --depth 1 https://git@bitbucket.org/genomicepidemiology/resfinder_db.git /resfinder_db && \
cd /resfinder_db && \
python3 INSTALL.py /usr/local/bin/kma_index non_interactive

# disinfinder database install
RUN mkdir /disinfinder_db && \
git clone --branch ${DISINFINDER_DB_COMMIT_HASH} --depth 1 https://git@bitbucket.org/genomicepidemiology/disinfinder_db.git /disinfinder_db && \
cd /disinfinder_db && \
python3 INSTALL.py /usr/local/bin/kma_index non_interactive

ENV PATH="${PATH}" \
# specifying database locations
CGE_RESFINDER_RESGENE_DB="/resfinder_db" \
CGE_RESFINDER_RESPOINT_DB="/pointfinder_db" \
CGE_DISINFINDER_DB="/disinfinder_db"

WORKDIR /data

