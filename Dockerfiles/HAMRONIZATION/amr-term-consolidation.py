# Term Harmonisation 
#When *hamronize* outputs the harmonized tsv, it does not harmonize the terms. So there are muliple rown that are the same gene at the same location identified by various tools and datasets.
#This program focuses on finding the equivalent genes and harmonizing the terms.
import pandas as pd
from collections import Counter
import sys
## Import Data
file = sys.argv[1]
hamr_output = pd.read_table(file)
print("PY file imported", file)


## Cleaning Data
#1. All the rows who's **analysis_software_name** is **resfinder**, and **'input_sequence_id'** is NA are collected in **read_resfinder**. Will be used later in XXX
#2. All the rows who's **analysis_software_name** is **resfinder**, and **'input_sequence_id'** is NOT are collected in **asm_resfinder**. The **'input_sequence_id'** will be manipulated to contain the true contig location of the gene. These are re-incorporated into **hamr_output**
#3. **hamr_output** is divided into hamr_idUNDER98 and hamr_idOVER98
def clean_df(hamr_output):
    # 1 NEW DF FROM RESFINDER THAT DOES NOT HAVE CONTIGS (READ). These do not have a input_sequence_id as it was generated from reads.
    read_resfinder = hamr_output[hamr_output['input_sequence_id'].isna()]
    
    # 2 NEW DF FROM RESFINDER THAT DOES HAVE CONTIGS (AMS). These do have a input_sequence_id because they come from contigs, but it is not in the correct format. 
    # eg "CCI165_S85_contig_8 length 181163 coverage 173.9 normalized_cov 0.95" becomes "CCI165_S85_contig_8"
    asm_resfinder = hamr_output[(hamr_output['analysis_software_name']=='resfinder') & (hamr_output['input_sequence_id'].notna())]
    asm_resfinder['input_sequence_id'] = asm_resfinder['input_sequence_id'].str.split().str[0]
    
    # All resfinder rows are removed from the dataframe  
    hamr_output = hamr_output[hamr_output['analysis_software_name']!='resfinder']
    
    # The modified resfinder rows from step 2 are added to the dataframe
    hamr_output = pd.concat([hamr_output,asm_resfinder])
    print("PY cleaned")
    return hamr_output, read_resfinder

# Generates hamr_output, where all input_sequence_id  are present and standarized; and res_resfinder which has empty input_sequence_id
# these will also be included into the harmonisation later on.
hamr_output, read_resfinder = clean_df(hamr_output)
# Function to create an inverse dictionary
def inverse_dictionary(data):
    inversed_data = {}
    for key, values in data.items():
        for value in values:
            inversed_data[value[0]] = key
    return inversed_data

# Function to harmonize terms based on optional cutoff values
def harmonize_terms(data, contig, cutoff_from=None, cutoff_to=None):
    harmonized_data = {}
    
    # If no cutoffs are provided
    if cutoff_from is None and cutoff_to is None:
        for key, values in data.items():
            value_collect = [value[0] for value in values]
            if value_collect:
                most_common_term = Counter(value_collect).most_common(1)[0][0]
                harmonized_data[most_common_term] = values

    # If both cutoff_from and cutoff_to are provided
    elif cutoff_from is not None and cutoff_to is not None:
        for key, values in data.items():
            value_collect = [value[0] for value in values if cutoff_from <= value[3] <= cutoff_to]
            if value_collect:
                most_common_term = Counter(value_collect).most_common(1)[0][0]
                harmonized_data[most_common_term] = values
    
    return harmonized_data

# Function to collect gene location information concisely
def concise_gene_location(data, row_collect, cutoff_from, cutoff_to):
    for key, value in data.items():
        count_in = 0
        count_all = len(value)
        
        score_all = [v[3] for v in value]
        score_in = [v[3] for v in value if cutoff_from <= v[3] <= cutoff_to]
        
        other_gene_symb_all = [v[0] for v in value]
        other_gene_symb_in = [v[0] for v in value if cutoff_from <= v[3] <= cutoff_to]

        avg_score_in = round(sum(score_in) / len(score_in),2) if score_in else None
        avg_score_all = round(sum(score_all) / len(score_all),2) if score_all else None

        count_in = len(score_in)
        
        row = [key, value[0][1], value[0][2], con, count_in, avg_score_in, other_gene_symb_in, count_all, avg_score_all, other_gene_symb_all]
        
        row_collect.append(row)
    
    return row_collect
# Collect unique contigs
contigs = hamr_output['input_sequence_id'].unique()
genome_inversed_harmonized_data_all = {}

row_collect = []
row_collect_all = []
for con in contigs:
    # Filter DataFrame for the current contig
    con_df = hamr_output[hamr_output['input_sequence_id'] == con]

    # Collect the following values
    gene_data = con_df[['gene_symbol', 'input_gene_start', 'input_gene_stop', 'sequence_identity']].values
    gene_index = con_df.index.values

    gene_collection = {}

    for gene, idx in zip(gene_data, gene_index):
        gene_sy, gene_sta, gene_sto, gene_id = gene
        range = 15
        start_1 = gene_sta - range
        start_2 = gene_sta + range
        stop_1 = gene_sta - range
        stop_2 = gene_sto + range

        key_found = False
        for k in gene_collection.keys():
            if (gene_sta >= k[0] and gene_sta <= k[1] and 
                gene_sto >= k[2] and gene_sto <= k[3]):
                gene_collection[k].append([gene_sy, gene_sta, gene_sto, gene_id, idx,con])
                key_found = True
                break
        
        if not key_found:
            gene_collection[(start_1, start_2, stop_1, stop_2)] = [[gene_sy, gene_sta, gene_sto, gene_id, idx,con]]

    # Process the gene collection for the current contig
    contig_harmonized_data_all = harmonize_terms(gene_collection, con)
    harmonized_data_over98 = harmonize_terms(gene_collection, con, cutoff_from= 98, cutoff_to=100)


    contig_inversed_harmonized_data_all = inverse_dictionary(contig_harmonized_data_all)

    row_collect = concise_gene_location(harmonized_data_over98, row_collect, 98,100)
    row_collect_all = concise_gene_location(contig_harmonized_data_all, row_collect_all, 98,100)
    
    genome_inversed_harmonized_data_all |= contig_inversed_harmonized_data_all

print("PY make harmonization files")
df_98 = pd.DataFrame(row_collect, columns=['harmonized_gene_symbol', 'start', 'end', "contig", "#hits>=98", "average_identity (>=98)", "Other_gene_symbols (>=98)","#hits (all)", "average_identity (all)", "Other_gene_symbols (all)"])
df_all = pd.DataFrame(row_collect_all, columns=['harmonized_gene_symbol', 'start', 'end', "contig", "#hits>=98", "average_identity (>=98)", "Other_gene_symbols (>=98)","#hits (all)", "average_identity (all)", "Other_gene_symbols (all)"])
hamr_output['harmonize_gene_symbol'] = hamr_output['gene_symbol'].map(genome_inversed_harmonized_data_all)
read_resfinder['harmonize_gene_symbol'] = read_resfinder['gene_symbol'].map(genome_inversed_harmonized_data_all)
hamr_output = pd.concat([hamr_output,read_resfinder])
hamr_output

# Filter: Sequence_Identity >= 98
hamr_IDisna = hamr_output[hamr_output['harmonize_gene_symbol'].isna()]
hamr_IDisna

# Assuming the `find_matches` function is defined as previously discussed
def find_matches(target, array):
    return [string for string in array if target[:4] in string]
    
harmonized_gene_symb = df_all['harmonized_gene_symbol'].tolist()

# Custom function to find matches for a row's third column value
def find_matches_for_row(row):
    return find_matches(row[2], harmonized_gene_symb)

# Apply the custom function to each row and store the result in a new column
hamr_IDisna['possible_harmonize_gene_symbol'] = hamr_IDisna.apply(find_matches_for_row, axis=1)



hamr_IDisna
hamr_IDisna.to_csv("hamronization_isna.tsv", sep='\t', index=False, mode='w')
hamr_output.to_csv("hamronization_all.tsv", sep='\t', index=False, mode='w')
df_98.to_csv("harmonized_amr_over98identity.tsv", sep = '\t', index = False, mode='w')
df_all.to_csv("harmonized_amr_allidentity.tsv", sep = '\t', index = False, mode='w')
print("PY DONE")
