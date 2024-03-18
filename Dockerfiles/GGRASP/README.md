# Welcome to the Dockerized GGRaSP

## Two ways of accessing the Docker

1. Building the image from the GGRaSP.Dockerfile, then running a container from the image 
2. Pulling the image from Dockerhub and running a container

## Note on running ggrasp

You can run ggrasp in two ways:
1. To allow for high-throughput analysis, an Rscript (ggrasp.R) file can run GGRaSP from the command line
2. On the R-console, ggrasp is available as a R-package.

## To test the Docker run the following script in a working directory where you want the outputs to be saved

### commandline method

mkdir /testrun && cd /testrun && Rscript /GGRaSP/ggrasp.R -i /GGRaSP/examples/Enter.ANI.mat -d 100 -o . --plotgmm -0 -1 -2 -3 -4 -5 -9 -7 -8

### R-console method

> library(ggrasp)
>enter.in.ggrasp <- ggrasp.load(system.file("extdata", "Enter.ANI.mat", package="ggrasp"), file.format="matrix", tree.method="complete", offset=100)
>enter.cluster <- ggrasp.cluster(enter.in.ggrasp);
>enter.cluster
>plot(enter.cluster, "gmm")

## Refer to

- https://www.jcvi.org/research/ggrasp
- https://github.com/JCVenterInstitute/ggrasp/
- https://cran.r-project.org/web/packages/ggrasp/index.html

GGRaSP: a R-package for selecting representative genomes using Gaussian mixture models. Thomas H Clarke, Lauren M Brinkac, Granger Sutton, and Derrick E Fouts. Bioinformatics, bty300, https://doi.org/10.1093/bioinformatics/bty300

