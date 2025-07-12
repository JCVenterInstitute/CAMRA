<p align="center">
  <img height="175" src="Images/camra-blue-logo.png">
</p>

# Combatting Antimicrobial Resistance in Africa Using Data Science

## Table of Contents
- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Docker Images](#docker-images)
- [WDL Workflow](#wdl-workflow)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## Introduction
### About CAMRA

CAMRA is a DS-I Africa research hub focused on **analyzing clinical and molecular data related to antimicrobial resistance (AMR)** in bacterial infections in Nigeria and Rwanda with the aim of **translating genomics of AMR to sensitive, rapid diagnostics, and effective therapeutics.**

### About CAMRA's Pipeline

A pipeline for the scientist who doesnt want to code. This pipeline is designed to assemble, qc, conduct extensive AMR detection with a viarety of tools, annotation and pangenome analysis, all to fulfill the scientific requirements for CAMRA and **with the intent to be used by other projects**.
The pipeline was designed to:
- run on Terra.bio (no commandline interaction, requires Terra.bio account)
- ran locally (requires familiarity with command line)
- Be modular. Users can decide to run entier pipeline. Or sections of the pipeline, example: only run the AMR worflow/section.
- Permit scientist with minimal computational skills the ability to process their bacterial isolates. 

## Prerequisites
### Terra.bio Run
- docker 
### Locall Run 



## What you will find here

The repository includes:

- Dockerfiles for JCVI-produced and 3rd party bioinformatic tools
- Workflows in the Workflow Description Language (WDL) tailored for Terra, a cloud-native research platform.

## Links

[CAMRA Website](https://camra.acegid.org/)

## Workflows

This pipeline is broken into four workflows: Assembly, Annotation, AMR Analysis, and Quality Control.

### Assembly

To run this workflow you need an active BV-BRC account.

#### Inputs

### Annotation

### AMR Analysis

Anti-microbial Resistance genes are identified in each assembly using multiple AMR programs along with the associated databases using docker images for each.
The results, despite the production of different output file formats, are harmonized with hAMRonize, an open source code that eases the downstream bioinformatic analysis of results. 

AMR Terms are consolidated and harmonized for each AMR Gene. 
RGI matching is prioritized in the reporting, but in instances where genes are only identified by another AMR, it will be used.


Tools used include:

- [amrfinder](https://hub.docker.com/r/ncbi/amr)
- [abricate](https://hub.docker.com/r/staphb/abricate)
- [resfinder](https://bitbucket.org/genomicepidemiology/resfinder/src/4.5.0/)
- [rgi](https://github.com/arpcard/rgi.git)
- [bvbrc](https://hub.docker.com/repository/docker/andrewrlapointe/bvbrc/general)
- [hamronization](https://hub.docker.com/r/finlaymaguire/hamronization)


### Quality Control

Tools used include:

- [MASH](https://mash.readthedocs.io/en/latest/)
- [MLST](https://github.com/tseemann/mlst)
- [entrezdirect](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [checkM](https://ecogenomics.github.io/CheckM/)
- [merqury](https://github.com/marbl/merqury)
- [quast](https://quast.sourceforge.net/)
- [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
