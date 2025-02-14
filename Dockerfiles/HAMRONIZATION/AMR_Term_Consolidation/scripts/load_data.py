import pandas as pd
import logging
import sys
from Bio import SeqIO
import obonet
import re


logger = logging.getLogger("term_consolidation")

# Step 1 - Load Data
def load_data(file) -> pd.DataFrame:
    """Load the data from the hamronize_amr_output.tsv."""
    try:
        hamr_output = pd.read_table(file)
        logger.info(f"> Total genes detected by all tools: {len(hamr_output)}")
        return hamr_output
    except FileNotFoundError:
        logger.error(f"File not found: {file}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading file {file}: {str(e)}")
        sys.exit(1)

def load_assembly(file) -> pd.DataFrame:
    """Load the data from the hamronize_amr_output.tsv."""
    try:
        assembly_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
        return assembly_dict
    except FileNotFoundError:
        logger.error(f"File not found: {file}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading file {file}: {str(e)}")
        sys.exit(1)

def retreive_card_ontology(file):
    """
    Load and process the CARD ontology graph from an OBO file.

    This function reads an ontology graph from the specified OBO file,
    extracts mappings of IDs to names, names to IDs, and synonyms to IDs,
    and cleans up the synonym terms for consistent use.

    Args:
        url (str): Path to the OBO file containing the CARD ontology.

    Returns:
        tuple: A tuple containing:
            - graph (networkx.MultiDiGraph): The ontology graph with nodes and edges.
            - id_to_name (dict): A dictionary mapping ontology IDs to their corresponding names.
            - name_to_id (dict): A dictionary mapping names to their corresponding ontology IDs.
            - synonym_to_id (dict): A dictionary mapping cleaned synonym terms to their corresponding ontology IDs.
    
    Notes:
        - Synonym terms are cleaned by removing metadata like 'EXACT' or 'CARD_Short_Name'.
        - All synonym terms are converted to lowercase for consistency.
        - The ontology graph is built using the `obonet` library.

    Example:
        >>> graph, id_to_name, name_to_id, synonym_to_id = retreive_card_ontology('/path/to/aro.obo')
        >>> print(id_to_name['ARO:0000001'])
        'Beta-lactamase'

    """
    # Load the ontology graph
    # url = '/Users/dmatute/Documents/CAMRA/lib/card_database/card_ontology/aro.obo'
    graph = obonet.read_obo(file)
    
    # Mapping of IDs to names
    id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
    
    # Mapping of names to IDs
    name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}
    
    # Dictionary to hold synonym to ID mapping
    synonym_to_id = {}
    
    # Iterate through the graph to clean synonyms and build the synonym_to_id mapping
    for node_id, data in graph.nodes(data=True):
        if 'synonym' in data:
            # Initialize the cleaned_synonyms list
            cleaned_synonyms = []
            
            # Clean each synonym in the list
            data['synonym'] = [re.sub(r'^"(.*?)" EXACT \[\]$', r'\1', synonym) for synonym in data['synonym']]
            
            # Clean any synonyms with 'EXACT' and 'CARD_Short_Name' as metadata
            data['synonym'] = [re.sub(r'^"(.*?)" EXACT CARD_Short_Name \[\]$', r'\1', synonym) for synonym in data['synonym']]
    
            # Lowercase all terms
            data['synonym'] = [synonym.lower() for synonym in data['synonym']]
            
            # Append cleaned synonyms to the list
            for synonym in data['synonym']:
                cleaned_synonyms.append(synonym)
            
            # Replace the original list with the cleaned list
            data['synonym'] = cleaned_synonyms
            
            # Map each cleaned synonym to the node_id
            for synonym in data['synonym']:
                synonym_to_id[synonym] = node_id

    return graph, id_to_name, name_to_id, synonym_to_id