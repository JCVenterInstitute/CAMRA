# Welcome to the Dockerized GGRaSP. 

## Two ways of accessing the Docker: 
1. Building the image from the GGRaSP.Dockerfile, then running a container from the image 
2. Pulling the image from Dockerhub and running a container

## To test the Docker run the following script in a working directory where you want the outputs to be saved: 
perl /LOCUST/locust/typer.pl -a /LOCUST/locust/example_data/Universal_Genes/alleles.fa --novel_schema -t raxml --biosample_list /LOCUST/locust/example_data/Universal_Genes/genomes_biosamples.txt --accession_list /LOCUST/locust/example_data/Universal_Genes/genomes_assemblies.txt --config /LOCUST/locust/default_config.ini

## Refer to: 
https://www.jcvi.org/research/ggrasp
https://github.com/JCVenterInstitute/ggrasp/
https://cran.r-project.org/web/packages/ggrasp/index.html

GGRaSP: a R-package for selecting representative genomes using Gaussian mixture models. Thomas H Clarke, Lauren M Brinkac, Granger Sutton, and Derrick E Fouts. Bioinformatics, bty300, https://doi.org/10.1093/bioinformatics/bty300

