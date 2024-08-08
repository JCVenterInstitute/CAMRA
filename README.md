<p align="center">
  <img height="175" src="Images/camra-blue-logo.png">
</p>

# Combatting Antimicrobial Resistance in Africa Using Data Science

## About

CAMRA is a DS-I Africa research hub focused on analyzing clinical and molecular data related to antimicrobial resistance (AMR) in bacterial infections in Nigeria and Rwanda with the aim of translating genomics of AMR to sensitive, rapid diagnostics, and effective therapeutics.

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

### Quality Control

Tools used include:

- [MASH](https://mash.readthedocs.io/en/latest/)
- [MLST](https://github.com/tseemann/mlst)
- [entrezdirect](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [checkM](https://ecogenomics.github.io/CheckM/)
- [merqury](https://github.com/marbl/merqury)
- [quast](https://quast.sourceforge.net/)
- [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
