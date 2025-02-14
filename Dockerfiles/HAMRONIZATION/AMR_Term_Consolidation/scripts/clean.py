from typing import Tuple
import pandas as pd
import logging

logger = logging.getLogger("term_consolidation")


# Step 2 - Clean Data


def clean_df(hamr_output: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Cleans and processes a hAMRonize output dataframe for consistent formatting and analysis.

    This function performs several cleaning and reformatting steps on the hAMRonize output:
    
    1. Extracts rows where `analysis_software_name` is `resfinder` and `input_sequence_id` is missing.
       - These rows are collected in `read_resfinder` and represent AMR hits generated from reads.
    2. Extracts rows where `analysis_software_name` is `resfinder` and `input_sequence_id` is present but improperly formatted.
       - These rows are collected in `asm_resfinder`, and the `input_sequence_id` column is reformatted to extract the true contig location.
       - The reformatted rows are added back to the dataframe.
    3. Extracts rows where `analysis_software_name` is `rgi`, and reformats the `input_sequence_id` column to extract the correct contig location.
       - Differentiates between `blast` and `diamond` hits based on the `input_file_name` column by updating `analysis_software_name` accordingly.
    4. Removes all original `resfinder` and `rgi` rows from the dataframe.
    5. Adds back the reformatted `asm_resfinder` and `rgi` rows.
    6. Removes duplicate rows for a clean final dataframe.

    Args:
        hamr_output (pd.DataFrame): The hAMRonize output dataframe containing AMR analysis results.
            Expected columns include:
            - `analysis_software_name`: Name of the analysis software used (e.g., 'resfinder', 'rgi').
            - `input_sequence_id`: Identifier for the sequence or contig.
            - `input_file_name`: Name of the input file used for the analysis.

    Returns:
        tuple:
            - pd.DataFrame: The cleaned and processed hAMRonize dataframe.
            - pd.DataFrame: A dataframe (`read_resfinder`) containing rows where `resfinder` was run on reads with no contig information.

    Notes:
        - Rows where `input_sequence_id` is reformatted include only the necessary identifiers for downstream analysis.
        - The function prints progress at each step to track changes in the dataframe size.
        - Duplicate rows are removed to ensure a clean dataset for further processing.

    Example:
        >>> cleaned_df, read_resfinder = clean_df(hamr_output)
        >>> print(cleaned_df.shape)
        (1000, 10)
        >>> print(read_resfinder.shape)
        (200, 10)
    """

    logger.info(
        f"CLEANING HAMRONIZE DATAFRAME ---------------------------------------------------------------------------------------------------\n")

    
    # 1. Extract rows where 'resfinder' analysis was run on reads (no contig information).
    # These rows have missing values in the 'input_sequence_id' column.
    # read_resfinder = hamr_output[hamr_output['input_sequence_id'].isna()]
    # logger.info(f"    Found {len(read_resfinder)} AMR hits for resfinder read analysis.")
    
    # 2. Extract rows where 'resfinder' analysis was run on assembled contigs.
    # These rows have improperly formatted 'input_sequence_id', which needs reformatting.
    # eg "CCI165_S85_contig_8 length 181163 coverage 173.9 normalized_cov 0.95" becomes "CCI165_S85_contig_8"
    asm_resfinder = hamr_output[
        (hamr_output['analysis_software_name']=='resfinder') & 
        (hamr_output['input_sequence_id'].notna())
        ]
    logger.info(f"    Found {len(asm_resfinder)} resfinder AMR hits that need contig name reformating.")

    # Reformat 'input_sequence_id' to extract the true contig name (e.g., "contig_8").
    asm_resfinder.loc[:, 'input_sequence_id'] = asm_resfinder['input_sequence_id'].str.split().str[0]
    

    # Step 3: Extract rows where 'rgi' analysis was run.
    # These rows also have improperly formatted 'input_sequence_id' that needs correction.
    # eg "CCI165_S85_contig_8_162"  becomes "CCI165_S85_contig_8"
    rgi = hamr_output[hamr_output['analysis_software_name'] == 'rgi']
    logger.info(f"    Found {len(rgi)} rgi AMR hits that need contig name reformating.")

    # Reformat 'input_sequence_id' to remove trailing identifiers (e.g., "_162").
    rgi.loc[:,'input_sequence_id'] = rgi['input_sequence_id'].str.split('_').str[:-1].apply(lambda x: '_'.join(x))


    # 4. Differentiate RGI hits by input file type ('blast' or 'diamond') for better categorization.
    # Update 'analysis_software_name' to indicate the type of hit.
    rgi.loc[rgi['input_file_name'].str.contains('blast', na=False), 'analysis_software_name'] = "rgi_blast"
    rgi.loc[rgi['input_file_name'].str.contains('diamond', na=False), 'analysis_software_name'] = "rgi_diamond"
    
    # 5. Remove all original 'resfinder' and 'rgi' rows from the dataframe
    # before re-adding the reformatted rows to ensure clean data.
    logger.info(f"    # Rows before deleting resfinder and rgi rows: {len(hamr_output)}")
    hamr_output = hamr_output[hamr_output['analysis_software_name']!='resfinder']
    hamr_output = hamr_output[hamr_output['analysis_software_name']!='rgi']
    logger.info(f"    # Rows after deleting resfinder and rgi rows: {len(hamr_output)}")
    
    # 6. Add the reformatted 'asm_resfinder' and 'rgi' rows back into the dataframe.
    hamr_output = pd.concat([hamr_output,asm_resfinder])
    hamr_output = pd.concat([hamr_output,rgi])
    logger.info(f"    # Rows after adding new input_sequence_names resfinder and rgi: {len(hamr_output)}")

    # 7. Remove any duplicate rows to ensure uniqueness.
    hamr_output.drop_duplicates(inplace=True)
    logger.info(f"    # Rows after removing duplicates: {len(hamr_output)}")

    logger.info("... Done cleaning.")

    # Return the cleaned dataframe and the extracted read-only resfinder rows.
    return hamr_output #, read_resfinder
