#!/bin/bash

# Define input file paths
HAMRONIZE_OUTPUT="/Users/dmatute/Projects/CAMRA/04_Scripts/CAMRA/Dockerfiles/HAMRONIZATION/AMR_Term_Consolidation/test_sample/hamronize_amr_output.tsv"
ONTOLOGY_FILE="/Users/dmatute/Projects/CAMRA/01_Raw_Data/External_Datasets/card_database/card_ontology/aro.obo"
ASSEMBLY_FILE="/Users/dmatute/Projects/CAMRA/04_Scripts/CAMRA/Dockerfiles/HAMRONIZATION/AMR_Term_Consolidation/test_sample/GCF_030128905.1/GCF_030128905.1_ASM3012890v1_genomic.fna"
DB_PROT_HOMOLOG="/Users/dmatute/Projects/CAMRA/01_Raw_Data/External_Datasets/card_database/card-data/makeblastdb_collection/protein_fasta_protein_homolog_model.tar.gz"
DB_PROT_VARIANT="/Users/dmatute/Projects/CAMRA/01_Raw_Data/External_Datasets/card_database/card-data/makeblastdb_collection/protein_fasta_protein_variant_model.tar.gz"
DB_NUCL_HOMOLOG="/Users/dmatute/Projects/CAMRA/01_Raw_Data/External_Datasets/card_database/card-data/makeblastdb_collection/nucleotide_fasta_protein_homolog_model.tar.gz"
DB_NUCL_VARIANT="/Users/dmatute/Projects/CAMRA/01_Raw_Data/External_Datasets/card_database/card-data/makeblastdb_collection/nucleotide_fasta_protein_variant_model.tar.gz"

# Run the Python script with the arguments
python /Users/dmatute/Projects/CAMRA/04_Scripts/CAMRA/Dockerfiles/HAMRONIZATION/AMR_Term_Consolidation/AMR_Term_Consolidation.py "$HAMRONIZE_OUTPUT" "$ONTOLOGY_FILE" "$ASSEMBLY_FILE" \
                    "$DB_PROT_HOMOLOG" "$DB_PROT_VARIANT" \
                    "$DB_NUCL_HOMOLOG" "$DB_NUCL_VARIANT"
