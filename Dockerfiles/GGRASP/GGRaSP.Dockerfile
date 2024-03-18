FROM --platform=linux/amd64 debian:rc-buggy
LABEL maintainer="Daniella Matute <dmatute@jcvi.org>" dockerfile-date-updated="Jan 11 2024" tool="GGRaSP" description="GGRaSP (Gaussian Genome Representative Selector with Prioritization), a R-package and associated executable Rscript program that generates a list of prioritized representative genomes from either supervised or unsupervised clustering of related genomes."

RUN apt-get update 
RUN apt-get install -y wget libssl-dev

#Install R
RUN apt-get autoclean && apt-get clean && apt-get update
RUN apt-get install -y dirmngr apt-transport-https ca-certificates software-properties-common gnupg2 
RUN apt-get install -y r-base 
RUN R -e "install.packages('getopt', repos = 'http://cran.us.r-project.org')"    
RUN R -e "install.packages('ape', repos = 'http://cran.us.r-project.org')"

# Installing ggrasp, the packages are in R. Do i need to make ggrasp.R an executable? 
RUN apt install -y cmake libcurl4-openssl-dev git
RUN R -e "install.packages(c('curl','httr','plotly','lme4','pbkrtest','car','bgmm'), dependencies = TRUE)"
RUN R -e "install.packages(c('mixtools','colorspace','methods'), dependencies = TRUE)"
RUN R -e "install.packages('ggplot2', version='0.9.1', dependencies = TRUE)"
RUN R -e "install.packages('ggrasp',dependencies = TRUE)"
RUN git clone https://github.com/JCVenterInstitute/GGRaSP.git


# Install MASH, it works, just type mash in the CL
RUN apt-get install -y mash -f

COPY README.md /README.md


