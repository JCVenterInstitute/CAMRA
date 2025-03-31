# Data handling and analysis
import pandas as pd
import numpy as np

# Utilities
import sys

# Modules
from scripts import card_matching,load_data,metadata, clean, utilities, term_consolidation
import logger_config

import numpy as np
import logging


# Main Script
def main(hamronize_output, ontology_file, assembly_file, db_paths):
    """Main execution function to orchestrate the workflow."""
    logger = logging.getLogger(__name__)

    # Step 1: Load data
    hamr_output = load_data.load_data(hamronize_output)
    assembly = load_data.load_assembly(assembly_file)
    graph, id_to_name, name_to_id, synonym_to_id = load_data.retreive_card_ontology(ontology_file)

    # Step 2: Clean data
    # TODO Figure out what to do with read_resfinder_df
    # hamr_output_df, read_resfinder_df = clean.clean_df(hamr_output)
    hamr_output_df = clean.clean_df(hamr_output)

    # Step 3: Group genes by contig
    hamr_output_df = utilities.group_genes(hamr_output_df)

    # Step 4: Generate metadata on tools and databases
    metadata.generate_metadata(hamr_output_df)
    print("#### db_paths", db_paths)
    # Step 5: Match AMR hits to CARD DB
    hamr_output_df = card_matching.match_to_card(
        hamr_output_df, 
        graph, 
        name_to_id, 
        synonym_to_id,
        id_to_name, 
        assembly,
        db_paths
        )
    
    matched_df = utilities.relevant_information(hamr_output_df)

    consolidated_terms_df =  term_consolidation.term_consolidation(matched_df, graph)
    # TODO term consolidation.
    
    # logger.info(
    #     f"... Success of CLAST\n {metadata.match_metadata(hamr_output_df,'card_match_type')}\n {hamr_output_df['card_match_type'].isna().sum()} np.nan AMR hits will be matched via BLASTp."
    # )
    
    logger.info(
        f"... Success of CLAST\n {metadata.match_metadata(matched_df,'card_match_type')}\n {matched_df['card_match_type'].isna().sum()} np.nan AMR hits will be matched via BLASTp."
    )

    logger.info("> Matching completed.")
    return matched_df, consolidated_terms_df  #hamr_output_df

if __name__ == "__main__":

    # LOGGING
    logger = logger_config.setup_logger()
    title = "AMR_Term_Consolidation started."
    logger.info(f"\n{'*' * (len(title) + 8)}\n* {title} *\n{'*' * (len(title) + 8)}")
    ###########################################################################

    # PARSING ARGUMENTS
    args = utilities.parse_arguments()
    # logger.info(f"> Arguments parsed: {args}")
    ###########################################################################

    # Store arguments in variables
    hamronize_file = args.hamronize_output_file
    logger.info(f"  > Hamronize file: {hamronize_file}")
    ontology_file = args.ontology_file
    logger.info(f"  > Ontology file: {ontology_file}")
    assembly_file = args.assembly_file
    logger.info(f"  > Assembly file: {assembly_file}")
    database_prot_homolog_file = args.database_prot_homolog_file
    logger.info(f"  > Database protein homolog file: {database_prot_homolog_file}")
    database_prot_variant_file = args.database_prot_variant_file
    logger.info(f"  > Database protein variant file: {database_prot_variant_file}")
    database_nucl_homolog_file = args.database_nucl_homolog_file
    logger.info(f"  > Database nucleotide homolog file: {database_nucl_homolog_file}")
    database_nucl_variant_file = args.database_nucl_variant_file
    logger.info(f"  > Database nucleotide variant file: {database_nucl_variant_file}")

    db_paths = {
    "prot_homolog": database_prot_homolog_file,
    "prot_variant": database_prot_variant_file ,
    "nucl_homolog": database_nucl_homolog_file,
    "nucl_variant": database_nucl_variant_file,
    }

    # logger.info(f"  > Database paths: {db_paths}")
    ###########################################################################

    # VALIDATE INPUT FILES
    logger.info(f"  > Validating input files...")
    utilities.validate_file(args.hamronize_output_file)
    utilities.validate_file(args.ontology_file)
    utilities.validate_file(args.assembly_file)
    # logger.info(f"> Input files validated: {args.hamronize_output_file}, {args.ontology_file}, {args.assembly_file}")

    for db_name, db_path in db_paths.items():
        utilities.validate_blast_db(db_path, db_name)
    ###########################################################################

    # RUN MAIN FUNCTION
    try:
        logger.info(f"  > Running main function...")
        hamronized_terms_df , consolidated_terms_df = main(
            hamronize_file, 
            ontology_file, 
            assembly_file, 
            db_paths
            )
        
        logger.info(f"PROCESSING COMPLETED SUCCESSFULLY ----------------------------------------------------------------------------------------------\n")
        logger.info(hamronized_terms_df.shape)
        logger.info(consolidated_terms_df)
        # Save the final dataframe as a TSV file
        hamronized_terms_df.to_csv("HARMONIZED_TERMS.tsv", sep="\t", index=False)
        consolidated_terms_df.to_csv("CONSOLIDATED_TERMS.tsv", sep="\t", index=False)
        logger.info(f"Processed data saved.")
        
    except Exception as e:
        logger.exception("An error occurred during processing.")
        sys.exit(1)

    logger.info("Script finished.")
