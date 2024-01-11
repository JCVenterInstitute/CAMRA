{\rtf1\ansi\ansicpg1252\cocoartf2709
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww24000\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Welcome to the Dockerized PanOct. \
\
## Two ways of accessing the Docker: \
1. Building the image from the PANOCT.Dockerfile, then running a container from the image \
2. Pulling the image from Dockerhub and running a container\
\
## To test the Docker run the following script in a working directory where you want the outputs to be saved: \
/panoct_v3.23/bin/panoct.pl -t /panoct_v3.23/example_dir/example_blast.txt -f /panoct_v3.23/example_dir/example_tags.txt -g /panoct_v3.23/example_dir/example.gene_att -P example_dir/example.pep -S Y -L 1 -M Y -H Y -V Y -N Y -F 1.33 -G y -c 0,25,50,75,100 -T \
\
## Refer to: \
https://www.jcvi.org/research/panoct\
https://pubmed.ncbi.nlm.nih.gov/22904089/\
https://sourceforge.net/p/panoct/home/Home/\
}