FROM --platform=linux/amd64 debian:rc-buggy
LABEL maintainer="Daniella Matute <dmatute@jcvi.org>" \
    dockerfile-date-updated="Mar 7 2024" \
    tool="entrez-direct" \
    description="Entrez Direct (EDirect) provides access to the NCBI's suite of interconnected databases (publication, sequence, structure, gene, variation, expression, etc.) from a Unix terminal window. "

RUN apt-get update 
RUN apt-get install -y wget cpanminus
RUN sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)" 
ENV PATH="/root/edirect:${PATH}"
