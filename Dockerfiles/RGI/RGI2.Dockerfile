FROM --platform=linux/amd64 continuumio/miniconda3

# Install system essentials
RUN apt-get update && apt-get install -y \
    wget \
    gcc && \
    apt-get clean && rm -rf /var/lib/apt/lists/* &&\
    mkdir /data

WORKDIR /opt

# Install RGI and dependencies
RUN wget https://github.com/arpcard/rgi/archive/refs/tags/6.0.3.tar.gz && \
    tar -xf 6.0.3.tar.gz && rm 6.0.3.tar.gz && \
    cd /opt/rgi-6.0.3 && \
    conda install -y \
        pip biopython=1.78 pytest pandas matplotlib seaborn \
        pyfaidx pyahocorasick pysam beautifulsoup4=4.9.3 requests \
        lxml blast=2.14.1 zlib prodigal=2.6.3 pyrodigal diamond=0.8.36 \
        oligoarrayaux=3.8 samtools bamtools=2.5.1 bedtools \
        bowtie2 bwa kma -c conda-forge -c bioconda && \
    pip install . 

    # Run RGI auto_load to initialize the environment
WORKDIR /opt/rgi-6.0.3
RUN rgi auto_load 

# Default working directory for data
WORKDIR /data