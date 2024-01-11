# Welcome to the Dockerized PanOct. 

## Two ways of accessing the Docker: 
1. Building the image from the PANOCT.Dockerfile, then running a container from the image 
2. Pulling the image from Dockerhub and running a container

## To test the Docker run the following script in a working directory where you want the outputs to be saved: 
/panoct_v3.23/bin/panoct.pl -t /panoct_v3.23/example_dir/example_blast.txt -f /panoct_v3.23/example_dir/example_tags.txt -g /panoct_v3.23/example_dir/example.gene_att -P example_dir/example.pep -S Y -L 1 -M Y -H Y -V Y -N Y -F 1.33 -G y -c 0,25,50,75,100 -T 

## Refer to: 
https://www.jcvi.org/research/panoct
https://pubmed.ncbi.nlm.nih.gov/22904089/
https://sourceforge.net/p/panoct/home/Home/


perl ../../typer.pl -a alleles.fa --novel_schema -t raxml --biosample_list genomes_biosamples.txt --accession_list genomes_assemblies.txt --config ../../default_config.ini
