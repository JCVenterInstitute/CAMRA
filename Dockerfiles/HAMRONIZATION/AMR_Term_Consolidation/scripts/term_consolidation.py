# Data handling and analysis
import pandas as pd
from typing import Tuple
from typing import List
from collections import Counter
import math
from scripts import load_data
import ast

class AMR_hit:
    def __init__(self, 
                 amr_tool_db,
                 amr_gene_symbol, 
                 amr_stop,
                 amr_start,
                 amr_contig,
                 amr_pident , 
                 amr_ref_accession,
                 match_name,
                 match_ref_accession,
                 match_pident):
        # Ensure amr_tool_db is treated as a tuple
        if isinstance(amr_tool_db, str):
            amr_tool_db = ast.literal_eval(amr_tool_db)
        def nan_to_str(field):
            if field:
                return field
            else:
                return ""

        self.amr_tool, self.amr_db = amr_tool_db  # Correct unpacking
        self.amr_gene_symbol = amr_gene_symbol
        self.amr_pident = amr_pident
        self.amr_ref_accession = amr_ref_accession
        self.amr_stop = amr_stop
        self.amr_start = amr_start
        self.amr_contig = amr_contig
        self.match_name = nan_to_str(match_name)
        self.match_ref_accession = nan_to_str(match_ref_accession)
        self.match_pident = nan_to_str(match_pident)
        
    def __repr__(self):
        return f"AMR_hit({self.amr_tool}, {self.amr_db},{self.amr_gene_symbol},{self.amr_pident},{self.amr_ref_accession},{self.amr_stop},{self.amr_start}, {self.amr_contig}, {self.match_name}, {self.match_ref_accession}, {self.match_pident})"

class loci_group:
    def __init__(self, group_id: int):
        self.group_id = group_id
        self.hits: List[AMR_hit] = []
        self.final_gene_info = None  # Store final results in a structured way

    def add_hit(self, hit: AMR_hit):
        self.hits.append(hit)

    def get_gene_info(self, graph, ref_acc):
        """Extract gene metadata from the graph."""
        if isinstance(ref_acc, str) and ref_acc in graph.nodes:
            node = graph.nodes[ref_acc]
            return {
                "namespace": node.get("namespace"),
                "definition": node.get("def"),
                "is_a": node.get("is_a"),
                "relationship": node.get("relationship"),
                "synonym": node.get("synonym"),
            }
        return {}

    def process_single_hit(self, graph):
        """Handles cases with a single hit."""
        hit = self.hits[0]
        self.final_gene_info = {
            "start": hit.amr_start,
            "stop": hit.amr_stop,
            "contig": hit.amr_contig,
            "hits": 1,
            "agreeing_hits": 1,
            "pident": hit.amr_pident if not math.isnan(hit.amr_pident) else 0,
            "ref_accession": hit.amr_ref_accession,
            "gene_name": hit.amr_gene_symbol,
            "database": hit.amr_db,
            "match_type": "Single Hit",
        }
        gene_info = self.get_gene_info(graph, hit.amr_ref_accession)

        if hasattr(hit, "match_ref_accession") and  (hit.match_ref_accession if isinstance(hit.match_ref_accession,str) else None) :
            match_metadata = self.parse_card_metadata(0, graph)
            gene_info.update(match_metadata)
            
        self.final_gene_info.update(gene_info)

    def process_multiple_hits(self, graph):
        """Handles cases with multiple hits by prioritizing 100% pident matches."""
        
        card_hits : Tuple[str, str] = [(hit.match_ref_accession, hit.match_pident) for hit in self.hits]
        card_hits = [(gene, float(pident) if isinstance(pident, str) and pident.strip() else 0) for gene, pident in card_hits]

        amr_hits : Tuple[str, str] = [(hit.amr_ref_accession, hit.amr_pident) for hit in self.hits]
        
        card_100 : List [str] = [ref for ref, pident in card_hits if pident == 100.0]
        amr_100 : List [str] = [ref for ref, pident in amr_hits if pident == 100.0]
        
        if card_100:
            self.select_most_common_hit(graph, card_100, "CARD", "Multiple hits")
        elif amr_100:
            self.select_most_common_hit(graph, amr_100, "AMR tool", "Multiple hits")
            
        else:
            card_highest_pident = max((pident if not math.isnan(pident) else 0) for _, pident in card_hits)
            # # card_highest_pident = max((float(pident) if not math.isnan(float(pident)) else 0) for _, pident in card_hits)
            # card_highest_pident = max(
            #     (float(pident) if pident.strip() and not math.isnan(float(pident)) else 0) 
            #     for _, pident in card_hits
            #     )

            card_highest_pident = card_highest_pident if not math.isnan(card_highest_pident) else 0
            card_high_identity_matches = [ref_acc for ref_acc, pident in card_hits if pident == card_highest_pident]


            amr_highest_pident = max((pident if not math.isnan(pident) else 0) for _, pident in amr_hits)
            amr_highest_pident = amr_highest_pident if not math.isnan(amr_highest_pident) else 0
            amr_high_identity_matches = [ref_acc for ref_acc, pident in amr_hits if pident == amr_highest_pident]
            
            if card_highest_pident  and card_high_identity_matches and (not amr_highest_pident or card_highest_pident >= amr_highest_pident)  :
                self.select_most_common_hit(graph, card_high_identity_matches, "CARD", "Multiple hits")

            elif amr_highest_pident and amr_highest_pident and (not card_highest_pident or amr_highest_pident >= card_highest_pident ):
                self.select_most_common_hit(graph, amr_high_identity_matches, "AMR tool", "Multiple hits")            
            
            
    def select_most_common_hit(self, graph, hits_list, db_name, consolidation_type):
        """Selects the most common reference accession from hits."""
        gene_counts = Counter(hits_list)
        
        if not gene_counts:  # Check if hits_list is empty
            print("Warning: No valid hits found.")
            self.final_gene_info = None
            return None
    
        most_common_gene, _ = gene_counts.most_common(1)[0]

        if db_name == "AMR tool":
            match_hit, card_count, most_common_index = self.find_card_match_from_amr(most_common_gene)
        else:
            match_hit, card_count, most_common_index = self.find_card_match_from_card(most_common_gene)
        
        if match_hit is None:
            print(f"Warning: No match found for {most_common_gene}. Returning None.")
            self.final_gene_info = None
            return None
            
        if db_name == "CARD":
            match_metadata = self.parse_card_metadata(most_common_index, graph) 
        elif db_name == "AMR tool":
            if self.hits[most_common_index].match_name:
                match_metadata = self.parse_card_metadata(most_common_index, graph)
            else :   
                match_metadata = {"pident":self.hits[most_common_index].amr_pident,
                                  "gene_name":self.hits[most_common_index].amr_gene_symbol,
                                  "database":self.hits[most_common_index].amr_db,
                                  "ref_accessiong":self.hits[most_common_index].amr_ref_accession}  
        self.final_gene_info = {
            "start": self.hits[most_common_index].amr_start,
            "stop": self.hits[most_common_index].amr_stop,
            "contig": self.hits[most_common_index].amr_contig,
            "match_type": consolidation_type,
            "hits" : len(self.hits),
            "agreeing_hits" : card_count
        }
        self.final_gene_info.update(match_metadata)
        return self.final_gene_info

    def find_card_match_from_amr(self, most_common_gene) -> tuple[str, int, int] | None:
        for index, record in enumerate(self.hits):
            if record.amr_ref_accession == most_common_gene:
                most_common_index = index       
        card_count = 0
        for i in self.hits:
            if i.match_ref_accession == self.hits[most_common_index].match_ref_accession:
                card_count += 1
        if card_count == 0:
            return None
        else:
            return self.hits[most_common_index].match_ref_accession , card_count , most_common_index

    def find_card_match_from_card(self, most_common_gene) -> tuple[str, int, int] | None:
        for index, record in enumerate(self.hits):

            try:
                match_pident = float(record.match_pident)
                if math.isnan(match_pident):
                    match_pident = 0  # Default to 0 if it's NaN
            except ValueError:
                match_pident = 0  # Default to 0 if conversion fails

            if record.match_ref_accession == most_common_gene and match_pident:
                most_common_index = index       
        card_count = 0
        for i in self.hits:
            if i.match_ref_accession == most_common_gene:
                card_count += 1
        if card_count == 0:
            return None
        else:
            return self.hits[most_common_index].match_ref_accession , card_count , most_common_index
            
    def parse_card_metadata (self, most_common_index, graph):
        ref_accession = self.hits[most_common_index].match_ref_accession

        try:
            match_pident = float(self.hits[most_common_index].match_pident)
            if math.isnan(match_pident):
                match_pident = 0  # Default to 0 if it's NaN
        except ValueError:
            match_pident = 0  # Default to 0 if conversion fails

        try:
            amr_pident = float(self.hits[most_common_index].amr_pident)
            if math.isnan(amr_pident):
                amr_pident = 0  # Default to 0 if it's NaN
        except ValueError:
            amr_pident = 0  # Default to 0 if conversion fails

        pident = match_pident if match_pident > amr_pident else amr_pident
        metadata_dict = {"pident":pident,
                         "gene_name":self.hits[most_common_index].match_name,
                         "database":"CARD",
                         "ref_accession":ref_accession}
        if ref_accession in graph.nodes:       
            graph_dict= {"name_space": graph.nodes[ref_accession]['namespace'] if 'namespace' in graph.nodes[ref_accession] and ref_accession and graph.nodes[ref_accession]['namespace'] else None,
                         "definition": graph.nodes[ref_accession]['def'] if 'def' in graph.nodes[ref_accession] and  ref_accession and graph.nodes[ref_accession]['def'] else None,
                         "is_a": graph.nodes[ref_accession]['is_a'] if 'is_a' in graph.nodes[ref_accession] and  ref_accession and graph.nodes[ref_accession]['is_a'] else None,
                         "relationship": graph.nodes[ref_accession]['relationship'] if  'relationship' in graph.nodes[ref_accession] and  ref_accession and graph.nodes[ref_accession]['relationship'] else None,
                         "synonym": graph.nodes[ref_accession]['synonym'] if 'synonym' in graph.nodes[ref_accession] and  ref_accession and graph.nodes[ref_accession]['synonym'] else None}
            metadata_dict.update(graph_dict) 
        return  metadata_dict
        
    def determine_final_ref_accession(self, graph):
        """Determines the final reference accession for the group."""
        if len(self.hits) == 1:
            self.process_single_hit(graph)
        else:
            self.process_multiple_hits(graph)

def create_loci_group_dataframe(loci_groups_dict, graph):
    data = []
    for group in loci_groups_dict:
        final_info = loci_groups_dict[group].final_gene_info
        if final_info is not None:  # Ensure it's not None before appending
            data.append(final_info)
    # Create a DataFrame from the list of dictionaries
    df = pd.DataFrame(data)
    return df


def term_consolidation (input_df, graph):
    # input_df = pd.read_csv("HARMONIZED_OUTPUT.tsv", sep = '\t')
    # input_df.head()

    loci_groups_dict = {}  # Dictionary to store LociGroup objects

    # Iterate each row in the dataframe 
    for _, row in input_df.iterrows():
        loci_id = row["loci_groups"]  # Loci group number
        
        # If the loci doesnt exsist in the dictionary, create it 
        if loci_id not in loci_groups_dict:
            loci_groups_dict[loci_id] = loci_group(loci_id)

        # Add the Amr hit to its respective loci in the dictionary
        hit = AMR_hit(
                    amr_tool_db = row["permutations"],
                    amr_gene_symbol = row["gene_symbol"],
                    amr_stop = row["input_gene_stop"],
                    amr_start = row["input_gene_start"],
                    amr_contig = row["input_sequence_id"],
                    amr_pident = row["sequence_identity"],
                    amr_ref_accession = row["reference_accession"],
                    match_name = row["card_match_name"],
                    match_ref_accession = row["card_match_id"],
                    match_pident = row["pident"],
        )

        loci_groups_dict[loci_id].add_hit(hit)

    for i in loci_groups_dict:
        loci_groups_dict[i].determine_final_ref_accession(graph)

    df = create_loci_group_dataframe(loci_groups_dict, graph)
    return df
