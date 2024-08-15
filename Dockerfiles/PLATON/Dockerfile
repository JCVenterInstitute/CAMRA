FROM debian:stable-slim


RUN apt update && \
apt upgrade -y


# Install CONDA
RUN apt -y install wget
RUN mkdir -p ~/miniconda3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
RUN bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
RUN rm -rf ~/miniconda3/miniconda.sh
RUN ~/miniconda3/bin/conda init bash

# Install Platon
RUN ~/miniconda3/bin/conda install -c conda-forge -c bioconda -c defaults platon

# Install DB
RUN wget https://zenodo.org/record/4066768/files/db.tar.gz
RUN tar -xzf db.tar.gz
RUN rm db.tar.gz

RUN chmod +x /root/miniconda3/bin/platon
ENV PATH=/root/miniconda3/bin/platon:$PATH


WORKDIR /data

CMD ["bash"]