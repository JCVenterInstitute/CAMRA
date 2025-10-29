
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
- run either:
	- on Terra.bio (no command line interaction, but requires Terra.bio account)
	- locally (requires familiarity with command line)
- Be modular. Users can decide to run entire pipeline or only sections of the pipeline, such as the AMR worflow/section.
- Permit scientist with minimal computational skills the ability to process their bacterial isolates. 

## Prerequisites
### Terra.bio Run
To run on [Terra.bio](https://app.terra.bio/#) requires both user account and an active project, both with can be setup according to the steps [HERE]
### Local Run 
- uses [Docker](https://www.docker.com/products/docker-desktop/) as a container for the different AMR programs
- can use [Cromwell](https://github.com/broadinstitute/cromwell/releases/tag/91)  from the Broad Institute to run terra commands

Description using exemplar data are shown [HERE]

## What you will find here

The repository includes:

- Dockerfiles for JCVI-produced and 3rd party bioinformatic tools
- Workflows in the Workflow Description Language (WDL) tailored for Terra, a cloud-native research platform. These can also be tested using Cromwell
- Example hAMRonization runs including:
   - (walk-throughs)[https://github.com/JCVenterInstitute/CAMRA/wiki/AMR-Analysis]
   - (input files)[https://github.com/JCVenterInstitute/CAMRA/tree/tclarke/examples/examples/Initial]
   - (output files)[https://github.com/JCVenterInstitute/CAMRA/tree/tclarke/examples/examples/Output]

## Links

[CAMRA Website](https://camra.acegid.org/)

## Workflows

This pipeline is broken into four workflows: Assembly, Annotation, AMR Analysis, and Quality Control.

### Assembly

To run this workflow you need an active BV-BRC account.

#### Inputs

[HERE]

### Annotation

Various Programs

### AMR Analysis

Anti-microbial Resistance genes are identified in each assembly using multiple AMR programs along with the associated databases using docker images for each. The programs used are:
- [amrfinder](https://hub.docker.com/r/ncbi/amr)
- [abricate](https://hub.docker.com/r/staphb/abricate)
- [resfinder](https://bitbucket.org/genomicepidemiology/resfinder/src/4.5.0/)
- [rgi](https://github.com/arpcard/rgi.git)
- [bvbrc](https://hub.docker.com/repository/docker/andrewrlapointe/bvbrc/general)

The results, despite the production of different output file formats, are harmonized with [hAMRonize](https://hub.docker.com/r/finlaymaguire/hamronization), an open source code that eases the downstream bioinformatic analysis of results.  RGI matching is prioritized in the reporting, but in instances where genes are only identified by another AMR, it will be used.

### Quality Control

Various tools are used as part of the pipeline to determine:

- the quality of the reads and the assembly
- taxonomic classification of the assembl


The programs are:
- [MASH](https://mash.readthedocs.io/en/latest/)
- [MLST](https://github.com/tseemann/mlst)
- [entrezdirect](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [checkM](https://ecogenomics.github.io/CheckM/)
- [merqury](https://github.com/marbl/merqury)
- [quast](https://quast.sourceforge.net/)
- [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
