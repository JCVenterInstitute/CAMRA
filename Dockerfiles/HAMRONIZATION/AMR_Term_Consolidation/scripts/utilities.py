import logging
from Bio.SeqRecord import SeqRecord
import pandas as pd
from Bio import SeqIO
import argparse
import sys
import os
import subprocess
import shutil

logger = logging.getLogger("term_consolidation")


def save_file_in_subdir(subdir, filename, content):
    """
    Saves a file in an intermediary directory inside the current working directory.

    Parameters:
    - subdir (str): Name of the subdirectory.
    - filename (str): Name of the file.
    - content (str): Text content to write.

    Returns:
    - str: Full path of the saved file.
    """
    # Get current working directory
    current_dir = os.getcwd()

    # Create the full path for the subdirectory
    dir_path = os.path.join(current_dir, subdir)

    # Ensure the directory exists
    os.makedirs(dir_path, exist_ok=True)

    # Create full file path
    file_path = os.path.join(dir_path, filename)

    # Write to the file
    with open(file_path, "w") as f:
        f.write(content)

    print(f"File saved at: {file_path}")
    return file_path

def validate_file(filepath):
    """Check if a file exists and is readable."""
    if not os.path.isfile(filepath):
        logger.error(f"    Missing file: {filepath}")
        sys.exit(1)
    else:
        logger.info(f"    Found file: {filepath}")

def validate_blast_db(db_path, db_name):
    def is_fasta_empty(fasta_file):
        """Check if a FASTA file is empty or contains no valid sequences."""
        try:
            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        return False  # Found a sequence header, so not empty
        except FileNotFoundError:
            print(f"Error: File '{fasta_file}' not found.")
            logger.error(f"    BLAST database is not found {db_path}")
            sys.exit(1)
        
        return True  # No valid sequences found

    def run_makeblastdb(fasta_file, db_name):
        """Run makeblastdb on the given FASTA file."""
        if 'nucl' in db_name:
            db_type = 'nucl'
        elif 'prot' in db_name:
            db_type = 'prot'
        cmd = [
            'makeblastdb',
            '-in', fasta_file,
            '-dbtype', db_type,
            '-out', db_name
        ]
        try:
            subprocess.run(cmd, check=True)
            logger.info(f"    BLAST database created successfuly: {db_name}")
            print(f"BLAST database created successfully: {db_name}")
        except subprocess.CalledProcessError as e:
            print(f"Error running makeblastdb: {e}")
            sys.exit(1)

    logger.info(f"    Checking BLAST database: {db_name} - {db_path}")

    if is_fasta_empty(db_path):
        print("FASTA file is empty or contains no valid sequences. Exiting.")
        logger.error(f"    BLAST database is empty: {db_path}")
        sys.exit(1)
    
    logger.info(f"    BLAST database is not empty: {db_path}")
    run_makeblastdb(db_path,db_name)
    """Check if a BLAST database is valid."""
    
    

    # # Ensure db_path is a valid string
    # if not isinstance(db_path, str) or not db_path:
    #     logger.error("    Invalid database path provided.")
    #     sys.exit(1)

    # # Expected BLAST database files based on the prefix
    # pin_file = db_path + ".pin" 
    # nin_file = db_path + ".nin"

    # if not os.path.isfile(pin_file) and not os.path.isfile(nin_file):
    #     logger.error(f"    BLAST database files missing for: {db_path}")
    #     sys.exit(1)

    # # Check if BLAST+ tools are installed
    # if not shutil.which("blastdbcmd"):
    #     logger.error("    BLAST tools not found. Ensure BLAST+ is installed and in PATH.")
    #     sys.exit(1)

    # # Validate BLAST database
    # try:
    #     result = subprocess.run(["blastdbcmd", "-db", db_path, "-info"], 
    #                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    #     if result.returncode != 0:
    #         logger.error(f"    Invalid BLAST database: {db_path}\n{result.stderr}")
    #         sys.exit(1)
    #     else:
    #         logger.info(f"    Validated BLAST database: {db_path}")
    # except Exception as e:
    #     logger.error(f"    Unexpected error during BLAST validation: {e}")
    #     sys.exit(1)
#     """Check if a BLAST database is valid."""
#     logger.info(f"    BLAST database: {db_path}")
#     logger.info(f"    BLAST database: {os.listdir(db_path)}")

#     # if not os.path.isfile(db_path + ".pin") and not os.path.isfile(db_path + ".nin"):
#     # if not os.path.isfile(os.path.join(db_path, db_path + ".pin")) and not os.path.isfile(os.path.join(db_path, db_path + ".nin")):
#     if not os.path.isfile(os.path.join(db_path, "db.pin")) and not os.path.isfile(os.path.join(db_path, "db.nin")):
    
#         logger.info(f"    > enter if")
#         logger.error(f"    BLAST database files missing for: {db_path}")
#         sys.exit(1)

#     try:
#         logger.info(f"    > enter try")
#         result = subprocess.run(["blastdbcmd", "-db", db_path, "-info"], 
#                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#         if result.returncode != 0:
#             logger.error(f"    Invalid BLAST database: {db_path}\n{result.stderr}")
#             sys.exit(1)
#         else:
#             logger.info(f"    Validated BLAST database: {db_path}")
#     except FileNotFoundError:
#         logger.info(f"    > enter except")
#         logger.error("    BLAST tools not found. Ensure BLAST+ is installed and in PATH.")
#         sys.exit(1)

def parse_arguments():
    """Parse command-line arguments for the AMR term consolidation script."""
    parser = argparse.ArgumentParser(description="Process AMR term consolidation.")

    # Define expected arguments
    parser.add_argument("hamronize_output_file", help="Path to the hamronize output file")
    parser.add_argument("ontology_file", help="Path to the ontology file")
    parser.add_argument("assembly_file", help="Path to the assembly FASTA file")
    parser.add_argument("database_prot_homolog_file", help="Path to protein homolog BLAST database")
    parser.add_argument("database_prot_variant_file", help="Path to protein variant BLAST database")
    parser.add_argument("database_nucl_homolog_file", help="Path to nucleotide homolog BLAST database")
    parser.add_argument("database_nucl_variant_file", help="Path to nucleotide variant BLAST database")

    return parser.parse_args()

# Step 3 - Group Loci
def group_genes (hamr_output_df):
    """
    Groups genes in the dataframe based on loci proximity and assigns group numbers.

    Parameters:
        hamr_output_df (pd.DataFrame): Input dataframe containing gene details.

    Returns:
        pd.DataFrame: Updated dataframe with a new column 'loci'.
    """

    logger.info(
        f"GROUPING DATAFRAME'S HITS BASED ON LOCI-----------------------------------------------------------------------------------------\n")
    
    # Extract working data values
    hamr_data = hamr_output_df[['input_sequence_id','gene_symbol', 
                                'input_gene_start', 'input_gene_stop']].values
    
    #list of dataframe indexes
    gene_index = hamr_output_df.index.values # Collect the indexes
    
    gene_groups = {}
    
    # Iterate through each of the AMR hits in the DataFrame
    for gene, idx in zip(hamr_data, gene_index):
        # Get each gene's start and stop loction and the contig. 
        gene_start, gene_stop = gene[2:4]
        gene_contig = gene[0]
    
        #Find the ranges of the start and stop
        range = 50
        
        start_1 = gene_start - range
        start_2 = gene_start + range
        stop_1 = gene_stop - range
        stop_2 = gene_stop + range
    
        match_found = False
        
        # Iterate the gene_group dictionary 
        for gene_key in gene_groups.keys():
            # If a loci match is found append the index to the dictionary. 
            if (gene_start >= gene_key[0] and gene_start <= gene_key[1] and 
                gene_stop >= gene_key[2] and gene_stop <= gene_key[3] and
                gene_contig == gene_key[4] ):
                gene_groups[gene_key].append(idx)
                match_found = True
                break
                
        # If a loci match is not found then create a new dictionayr key and let the value be the index.    
        if not match_found:
            gene_groups[(start_1, start_2, stop_1, stop_2, gene_contig)] = [idx]
    
    # Change keys to numbers from 0 to len(dictionary) - 1
    gene_groups = {i: v for i, (k, v) in enumerate(gene_groups.items())}
    
    # Iterate through the dictionary and assign group numbers to the dataframe
    for group_number, indices in gene_groups.items():
        hamr_output_df.loc[indices, 'loci_groups'] = group_number
    hamr_output_df['loci_groups'] = hamr_output_df['loci_groups'].astype(int)
    logger.info(f"... Done grouping {len(hamr_output_df)} AMR hits of various tools into {len(gene_groups)} loci groups.")
    return hamr_output_df

# Step 6 - BLASTp unmached hits. 
def extract_sequence(assembly:dict, contig:str, contig_header:str, start:int, end:int, gene_name:str, count : int) :
    if start < end:
        nuc_sequence = assembly[contig_header].seq[start:end]
        pro_sequence = nuc_sequence.translate(table=11,to_stop=True)
        pro_record = SeqRecord(seq=pro_sequence, id=gene_name, description=f"{contig}:{start}-{end}")
        nuc_record = SeqRecord(seq=nuc_sequence, id=gene_name, description=f"{contig}:{start}-{end}")

    elif start > end:
        nuc_sequence = assembly[contig_header].seq[end:start]
        nuc_sequence = nuc_sequence.reverse_complement()
        pro_sequence = nuc_sequence.translate(table=11)
        pro_record = SeqRecord(seq=pro_sequence, id=gene_name, description=f"{contig}:{end}-{start}")
        nuc_record = SeqRecord(seq=nuc_sequence, id=gene_name, description=f"{contig}:{end}-{start}")

    if pro_record and nuc_record:
        logger.info(f"    {count}.successful seq extract : {contig} {start}:{end}")
        count = count + 1
        return pro_record, nuc_record, count
    elif pro_record or nuc_record: 
        logger.info(f"    {count}.semi-successful seq extract : {contig} {start}:{end}")
        count = count + 1
        return pro_record, nuc_record, count
    else:
        logger.info(f"    {count}.failed seq extract : {contig} {start}:{end}")
        count = count + 1
        return None, None, None

def update_df(df, subset_df, match_type) -> pd.DataFrame:

    # Ensure 'pident' column exists in df
    if 'pident' not in df.columns:
        df['pident'] = float()# Initialize with very low values
    
    # Iterate over rows in subset_df
    for idx, row in subset_df.iterrows():
        if idx in df.index:
            # Compare 'pident' values
            if (
                row['pident'] > df.at[idx, 'pident'] and 
                str(df.at[idx, 'card_match_type']) not in ['gene_symbol(SYNONYM)', 'gene_name(SYNONYM)', 'reference_accession', 'gene_symbol', 'gene_name']
            ):
                df.at[idx, 'card_match_id'] = row['reference_number']
                df.at[idx, 'card_match_name'] = row['reference_name']  # Fixed typo
                df.at[idx, 'pident'] = row['pident']
                df.at[idx, 'card_match_type'] = match_type

        else:
            # If index does not exist in df, add it
            df.loc[idx] = row[['reference_number', 'reference_name', 'pident']]
            df.at[idx, 'card_match_type'] = match_type
    
    return df


def make_fasta_file(df: pd.DataFrame, assembly: dict, output_dir: str, output_prefix: str) -> None:
    """
    Extracts AMR sequences from an assembly and saves them as FASTA and TXT files.

    Parameters:
    - df (pd.DataFrame): DataFrame containing AMR gene information.
    - assembly (dict): Dictionary containing contig sequences.
    - output_dir (str): Directory to save output files.
    - output_prefix (str): Prefix for the output files.

    Returns:
    - None (Saves files in the specified directory)
    """
    
    count = 1
    logger.info(
        f"EXTRACTING AMR SEQUENCES FROM ASSEMBLY------------------------------------------------------------------------------------------\n")

    nucleotide_records = []
    protein_records = []
    metadata_list = []

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    for _, row in df.iterrows():
        contig = row["input_sequence_id"]
        start = row["input_gene_start"]
        end = row["input_gene_stop"]
        gene_name = f"{row.name}_{'_'.join(row['permutations'])}"

        # Extract the contig sequence from assembly dictionary
        for contig_header in assembly.keys():
            if contig in contig_header:
                pro_record, nuc_record, count = extract_sequence(assembly, contig, contig_header, start, end, gene_name, count)
                nucleotide_records.append(nuc_record)
                protein_records.append(pro_record)

                # Collect metadata
                metadata_list.append(f"{gene_name}\t{contig}\t{start}\t{end}")
                break

    # Define file paths
    nucleotide_fasta_path = os.path.join(output_dir, f"{output_prefix}_nucleotide.fasta")
    protein_fasta_path = os.path.join(output_dir, f"{output_prefix}_protein.fasta")
    metadata_txt_path = os.path.join(output_dir, f"{output_prefix}_metadata.txt")

    # Save nucleotide and protein FASTA files
    if nucleotide_records:
        SeqIO.write(nucleotide_records, nucleotide_fasta_path, "fasta")
        logger.info(f"  > Nucleotide FASTA saved at {nucleotide_fasta_path}")

    if protein_records:
        SeqIO.write(protein_records, protein_fasta_path, "fasta")
        logger.info(f"  > Protein FASTA saved at {protein_fasta_path}")

    # Save metadata TXT file
    if metadata_list:
        with open(metadata_txt_path, "w") as f:
            f.write("Gene_ID\tContig\tStart\tEnd\n")  # Add header
            f.write("\n".join(metadata_list))
        logger.info(f"  > Metadata TXT saved at {metadata_txt_path}")

def save_intermediary_file(output_dir: str, filename: str, data):
    """
    Saves a Pandas DataFrame or raw BLAST output into an intermediary file.

    Parameters:
    - output_dir (str): Directory where the file will be saved.
    - filename (str): Name of the output file (with extension).
    - data (pd.DataFrame or str): Data to be saved. Can be a Pandas DataFrame or raw text.

    Returns:
    - str: Full path of the saved file.
    """

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define full file path
    file_path = os.path.join(output_dir, filename)

    # Handle saving based on data type
    if isinstance(data, pd.DataFrame):
        # Save DataFrame in appropriate format
        if filename.endswith(".csv"):
            data.to_csv(file_path, index=False)
        elif filename.endswith(".tsv"):
            data.to_csv(file_path, index=False, sep="\t")
        elif filename.endswith(".xlsx"):
            data.to_excel(file_path, index=False, engine="openpyxl")
        else:
            raise ValueError("Unsupported file format. Use .csv, .tsv, or .xlsx for DataFrames.")
    elif isinstance(data, str):
        # Save raw BLAST output as a text file
        with open(file_path, "w") as f:
            f.write(data)
    else:
        raise TypeError("Data must be a Pandas DataFrame or a string.")

    print(f"File saved at: {file_path}")
    return file_path


def relevant_information(hamr_output_df):
    
    missing_cols = [col for col in ['loci_groups', 'permutations', 'gene_symbol', 'sequence_identity', 
                                'reference_accession', 'card_match_type', 'card_match_name', 
                                'card_match_id', 'pident'] if col not in hamr_output_df.columns]

    if missing_cols:
        print(f"Warning: The following columns are missing: {missing_cols}")
    else:
        df = hamr_output_df[['loci_groups', 'permutations', 'gene_symbol', 'sequence_identity', 
                                        'reference_accession', 'card_match_type', 'card_match_name', 
                                        'card_match_id', 'pident','input_gene_start','input_gene_stop','input_sequence_id']]
        df["pident"] = df["pident"].replace(0.0, "")
        return df