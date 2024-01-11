FROM --platform=linux/amd64 debian:rc-buggy
LABEL maintainer="Daniella Matute <dmatute@jcvi.org>"

RUN apt-get update 
RUN apt-get install -y wget libssl-dev

#Install R
RUN apt-get autoclean && apt-get clean && apt-get update
RUN apt-get install -y dirmngr apt-transport-https ca-certificates software-properties-common gnupg2 
RUN apt-get install -y r-base 
RUN R -e "install.packages('getopt', repos = 'http://cran.us.r-project.org')"    
RUN R -e "install.packages('ape', repos = 'http://cran.us.r-project.org')"

# Installing ggrasp, the packages are in R. Do i need to make ggrasp.R an executable? 
RUN apt install -y cmake libcurl4-openssl-dev
RUN R -e "install.packages(c('curl','httr','plotly','lme4','pbkrtest','car','bgmm'), dependencies = TRUE)"
RUN R -e "install.packages(c('mixtools','colorspace','methods'), dependencies = TRUE)"
RUN R -e "install.packages('ggplot2', version='0.9.1', dependencies = TRUE)"
RUN R -e "install.packages('ggrasp',dependencies = TRUE)"
COPY GGRaSP-master.zip /
RUN unzip GGRaSP-master.zip && rm -rd GGRaSP-master.zip


# Install MASH, it works, just type mash in the CL
RUN apt-get install -y mash -f


