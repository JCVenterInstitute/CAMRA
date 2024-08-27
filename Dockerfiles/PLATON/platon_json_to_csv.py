import json
import pandas as pd
import sys

# Function that checks what AMR hits exsist in the plasmid hits
def amr_in_plasmid(plasmid_df, amr_df) :
    all_contigs = plasmid_df['contig_id'].unique()
    plasmid_amr = []
    
    
    #Iterate through contigs that have plasmid hits
    for contig in all_contigs:
        plasmid_info = plasmid_df.loc[plasmid_df['contig_id'] == contig]
        start, end = None, None
        # This assumes that each contig can only match with a single plasmid
        # Define where in the contig is the plasmid hit
        start = plasmid_info['contig_start'].iloc[0]
        end = plasmid_info['contig_end'].iloc[0]
        plasmid_id = plasmid_info['plasmid_id'].iloc[0]
        # start = list(plasmid_df[plasmid_df['contig_id']== contig]['contig_start'])[0]
        # end = list(plasmid_df[plasmid_df['contig_id']== contig]['contig_end'])[0]
        # plasmid_id = list(plasmid_df[plasmid_df['contig_id']== contig]['plasmid_id'])[0]
        if start is not None and end is not None:
            # Dataframe collects AMR hits that are in Plamids
            df = amr_df[(amr_df['contig_id'] == contig) & 
                            (amr_df['start'] > start) & 
                            (amr_df['start'] < end) & 
                            (amr_df['end'] > start) & 
                            (amr_df['end'] < end)]
            df['plasmid_id'] = plasmid_id
            
            # Append the filtered AMR hits to the list
            plasmid_amr.append(df)
    
    # Concatenate all filtered DataFrames into a single DataFrame
    plasmid_amr_df = pd.concat(plasmid_amr)
    return plasmid_amr_df
                    
# Extracts AMR, Plasmid and ORFs from json file produced by Palton 1.7
def extract_json(json_file, mode):
    plasmid_hits = []
    amr_hits = []
    orfs_hits = []
    
    # Opening JSON file
    f = open(json_file)
    
    # returns JSON object as 
    # a dictionary
    genome = json.load(f)
    
    # Iterating through the json
    # list
    for contig in genome:
        if genome[contig]['plasmid_hits']:
            for hit_i in range(len(genome[contig]['plasmid_hits'])):
                hit = genome[contig]['plasmid_hits'][hit_i]
                hit['plasmid_id'] = hit['plasmid']['id']
                hit['plasmid_length'] = hit['plasmid']['length']
                hit['plasmid_hit_size'] = abs(hit['plasmid_end'] - hit['plasmid_start'])
                hit['contig_hit_size'] = abs(hit['contig_end'] - hit['contig_start'])
                hit['plasmid_percentage']= hit['plasmid_hit_size']/hit['plasmid_length']*100
                hit['contig_id'] = contig
                hit['contig_length']= genome[contig]['length']
                hit['mode']= mode
                hit['contig_hit_percentage'] = hit['contig_hit_size']/hit['contig_length']*100
                del hit['plasmid']
                plasmid_hits.append(hit)

        if genome[contig]['amr_hits']:
            for hit_i in range(len(genome[contig]['amr_hits'])):
                hit = genome[contig]['amr_hits'][hit_i]
                hit['contig_id'] = contig
                hit['contig_length']= genome[contig]['length']
                hit['mode']= mode
                amr_hits.append(hit)
                
        if genome[contig]['orfs']:
            for hit_i in genome[contig]['orfs']:
                hit = genome[contig]['orfs'][hit_i]
                hit['contig_id'] = contig
                hit['contig_length']= genome[contig]['length']
                hit['mode']= mode
                orfs_hits.append(hit)
    # Closing file
    f.close()

    # Make Dataframes
    plasmid_df = pd.DataFrame(plasmid_hits)
    amr_df = pd.DataFrame(amr_hits)
    orfs_df = pd.DataFrame(orfs_hits)

    # Rearange Columns
    plasmid_df = plasmid_df.loc[:,['contig_id','contig_length', 'plasmid_id','plasmid_length','coverage', 'identity','contig_hit_size', 'plasmid_hit_size' ,'contig_hit_percentage','plasmid_percentage','contig_start', 'contig_end', 'plasmid_start', 'plasmid_end','mode']]

    # Call function that checks what AMR hits exsist in the plasmid hits
    plasmid_amr = amr_in_plasmid(plasmid_df, amr_df)
    
    return plasmid_df, amr_df, orfs_df, plasmid_amr
    

file = sys.argv[1]

# Platon Mode: Accuracy 
genome_plasmid, genome_amr, genome_orfs, plasmid_amr_df = extract_json(file, "accuracy")


genome_plasmid.to_csv('./Platon_Output/plasmid_hits.csv', index=False, mode='w')  
genome_amr.to_csv('./Platon_Output/amr_hits.csv', index=False, mode='w')  
genome_orfs.to_csv('./Platon_Output/orfs_hits.csv', index=False, mode='w')  
plasmid_amr_df.to_csv('./Platon_Output/amr_in_plasmid.csv', index=False, mode='w')  