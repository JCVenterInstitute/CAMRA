import logging
import pandas as pd
import numpy as np
import subprocess
from scripts import metadata, utilities
import os


logger = logging.getLogger("term_consolidation")


def match_to_card(hamr_output_df, graph, name_to_id, synonym_to_id, id_to_name, assembly, db_paths):
    """
    Matches AMR hits to CARD DB and attempts BLAST-based matching for unmatched hits.
    
    Parameters:
    - hamr_output_df (pd.DataFrame): DataFrame with AMR gene matches.
    - graph, name_to_id, synonym_to_id, id_to_name: CARD ontology mapping.
    - assembly (dict): Genome assembly sequences.
    - db_paths (dict): Dictionary containing paths to BLAST databases.
    
    Returns:
    - pd.DataFrame: Updated DataFrame with matched AMR genes.
    """

    logger.info(
        f"MATCHING AMR HITS TO CARD DB----------------------------------------------------------------------------------------------------------\n")

    output = f"intermediary"
    # Ensure output directory exists
    os.makedirs(output, exist_ok=True)

    # Initialize match columns
    hamr_output_df[["card_match_type", "card_match_name", "card_match_id"]] = np.nan

    # Apply CARD matching function
    hamr_output_df = hamr_output_df.apply(
        matching_to_card, 
        args=(graph, name_to_id, synonym_to_id, id_to_name), 
        axis=1
    )

    # Log match success
    logger.info(
        f"... Success of Matches:\n {metadata.match_metadata(hamr_output_df,'card_match_type')}\n"
        f"   {hamr_output_df['card_match_type'].isna().sum()} AMR hits np.nan hits did not surpass the pident >= 75."
    )

    # Collect and generate metadata for missing matches
    missing_matches_df = hamr_output_df[hamr_output_df['card_match_id'].isna()]
    if not missing_matches_df.empty:
        missing_matches_df.to_csv(f"{output}/missing_matches.tsv", sep="\t", index=True)
        metadata.generate_metadata(missing_matches_df)
        utilities.make_fasta_file(missing_matches_df, assembly, "intermediary","missing_matches")

        # Define BLAST parameters
        blast_tasks = [
            ("./intermediary/missing_matches_protein.fasta", "blastp", "prot_homolog", "card_blastp_homolog"),
            ("./intermediary/missing_matches_protein.fasta", "blastp", "prot_variant", "card_blastp_variant"),
            ("./intermediary/missing_matches_nucleotide.fasta", "blastx", "prot_homolog", "card_blastx_homolog"),
            ("./intermediary/missing_matches_nucleotide.fasta", "blastx", "prot_variant", "card_blastx_variant"),
            ("./intermediary/missing_matches_nucleotide.fasta", "blastn", "nucl_homolog", "card_blastn_homolog"),
            ("./intermediary/missing_matches_nucleotide.fasta", "blastn", "nucl_variant", "card_blastn_variant"),
        ]

        logger.info(f"BLASTING MISSING HITS ----------------------------------------------------------------------------------------------------------\n")
        # Run BLAST for each task and update dataframe
        for fasta_file, blast_type, db_path, match_type in blast_tasks:
            print("####BLAST_TASKS", fasta_file, blast_type, db_path, match_type)
            try:
                blast_results = blast_missing_hits(fasta_file, db_path, match_type, blast_type)
                hamr_output_df = utilities.update_df(hamr_output_df, blast_results, match_type)
            except Exception as e:
                logger.error(f"BLAST {blast_type} failed for {match_type}: {e}")

        logger.info(
        f"... Success of np.nan Matches:\n {metadata.match_metadata(hamr_output_df,'card_match_type')}\n"
        f"   {hamr_output_df['card_match_type'].isna().sum()} np.nan AMR hits will be matched via BLASTp."
    )

    return hamr_output_df

def matching_to_card(row, graph, name_to_id, synonym_to_id, id_to_name):
    """
    Matches a gene symbol to its corresponding CARD identifier using multiple matching strategies.

    This function checks for exact matches and matches based on gene symbols, gene_name, synonyms, 
    and reference accessions. The function also ensures that matches are 
    added to the `row` with relevant match details.

    Parameters:
    ----------
    row : pandas.Series
        A row from a dataframe containing the following columns:
        - `reference_accession`: The reference accession for the gene.
        - `gene_symbol`: The gene symbol to be matched.
        - `gene_name`: The gene name to be matched.
    graph : dict
        A dictionary containing reference accessions as keys and their corresponding CARD identifiers as values.
    name_to_id : dict
        A dictionary mapping gene names to their corresponding CARD identifiers.
    synonym_to_id : dict
        A dictionary mapping gene synonyms to their corresponding CARD identifiers.
    id_to_name : dict
        A dictionary mapping CARD identifiers to their corresponding gene_names.

    Returns:
    -------
    pandas.Series
        The updated `row` with the following additional columns:
        - `card_match_name`: The matched gene symbol or accession name.
        - `card_match_id`: The CARD identifier corresponding to the match.
        - `card_match_type`: The type of match ('reference_accession', 'gene_symbol', 'gene_symbol(SYNONYM)', or 'requires_blast').

    Notes:
    -----
    - Exact matches are prioritized.

    """
    match_found = False # Flag to track if a match has been found

    # 1. Check for exact match to reference_accessionn
    if "ARO:" + row['reference_accession'] in graph and not match_found:
        row["card_match_name"] = id_to_name["ARO:" + row['reference_accession']]
        row["card_match_id"] = "ARO:" + row['reference_accession']
        row["card_match_type"]  = 'reference_accession'
        match_found = True

    # 2. Match gene_symbol to name
    if row['gene_symbol'].lower() in name_to_id and not match_found:
        row["card_match_name"] = row['gene_symbol']
        row["card_match_id"] = name_to_id[row['gene_symbol']]
        row["card_match_type"] = 'gene_symbol'
        match_found = True
        
    # 3. Match gene_symbol to synonym
    if row['gene_symbol'].lower() in synonym_to_id and not match_found:
        row["card_match_name"] = row['gene_symbol']
        row["card_match_id"] = synonym_to_id[row['gene_symbol'].lower()]
        row["card_match_type"] = 'gene_symbol(SYNONYM)'
        match_found = True
    
    # 4. Match gene_symbol to name
    if row['gene_name'].lower() in name_to_id and not match_found:
        row["card_match_name"] = row['gene_name']
        row["card_match_id"] = name_to_id[row['gene_name']]
        row["card_match_type"] = 'gene_name'
        match_found = True

    # 5. Match gene_symbol to synonym
    if row['gene_name'].lower() in synonym_to_id and not match_found:
        row["card_match_name"] = row['gene_name']
        row["card_match_id"] = synonym_to_id[row['gene_name'].lower()]
        row["card_match_type"] = 'gene_name(SYNONYM)'
        match_found = True   
        
    return row     

def blast_missing_hits(fasta_file:str, database: str, output_file : str, blast_type:str) -> pd.DataFrame:
    num_reads = sum(1 for line in open(fasta_file) if line.startswith(">"))
    logger.info(
        f"> {blast_type.upper()} {num_reads} AMR hits against {output_file} DB")

    # Run BLAST with additional filtering options
    blast_command = [
        blast_type,  # Change to "blastp" for protein sequences
        "-query", fasta_file,
        "-db", database,
        "-out", f"intermediary/{output_file}_results.txt",
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-num_threads", "8",
        "-max_hsps", "1",  # Ensure only one high-scoring segment per subject
    ]
    
    subprocess.run(blast_command, check=True)

    # Read BLAST output
    blast_df = pd.read_csv(f"intermediary/{output_file}_results.txt", sep="\t", header=None, names=[
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ])
    # Extract the number before the first "_" and set it as the index
    blast_df.index = blast_df["qseqid"].str.extract(r"^(\d+)")[0].astype(int)
    # Filter for hits with pident > 75
    blast_df = blast_df[blast_df["pident"] > 75]
    blast_df = blast_df.sort_values(by=['qseqid'],ascending=False)

    # Ensure only one hit per query sequence (qseqid)
    blast_df = blast_df.drop_duplicates(subset="qseqid", keep="first")

    # Extract the reference number (digits after 'ARO:')
    # blast_df['reference_number'] = blast_df['sseqid'].str.extract(r'ARO:(\d+)')
    # blast_df['reference_number'] = blast_df['sseqid'].str.split("|")[4]
    blast_df['reference_number'] = blast_df['sseqid'].str.extract(r'(ARO:\d+)')


    # Extract the reference name (text after last '|')
    blast_df['reference_name'] = blast_df['sseqid'].str.extract(r'\|([^|]+)$')
    logger.info("    ... BLAST Complete.")

    
    return blast_df

