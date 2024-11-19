import pandas as pd
import json
import re
import numpy as np
import sys

# READ FILES
quality_json_file = sys.argv[1]
genome_annotation_file = sys.argv[2]
amr_json = sys.argv[3]

with open(genome_annotation_file, 'r') as file:
        annotation_genome = json.load(file)

with open(quality_json_file, 'r') as file:
    quality_json = json.load(file)

# CREATE QUALITY DATAFRAME - THIS DATAFRAME CONTAINS THE LIST OF AMR WITH HIGH CONFIDENCE, BUT IT CONTAINS LITTLE ANNOTATION INFORMATION
quality_json_df = pd.DataFrame(quality_json['amr_genes'], columns=['feature_id', 'gene_symbol', 'gene_definition', 'function'])

# CREATE FEATURE DATAFRAME - HAS ALL THE FEATURES AND DETAILS OF THE ANNOTATION EVENT, BUT NOT THE ANNOTATION TOOL
# Make dataframe with AMR features
features_dict = {}
# Define the pattern where * can be A-Z or 0-9. We will use this to find annotation_event_ids
pattern_id = r'^[A-Z0-9]{8}-[A-Z0-9]{4}-[A-Z0-9]{4}-[A-Z0-9]{4}-[A-Z0-9]{12}$'
# feature number
count = 0

# Annotate through all the annoated features in annotation.genome
for i in annotation_genome['features']:
    ID = i['id']
    # If the annotated feature id is in the quality.json['feature_id'] column, then add it to the df
    if ID in quality_json_df['feature_id'].tolist(): 
        features_dict[ID]=i
        features_dict[ID]['feature_number'] = count
        # FIX COLUMNS : location
        features_dict[ID]['contig'] = features_dict[ID]['location'][0][0]
        features_dict[ID]['start'] = features_dict[ID]['location'][0][1]
        features_dict[ID]['orientation'] = features_dict[ID]['location'][0][2]
        features_dict[ID]['end'] = features_dict[ID]['location'][0][3]
        # FIX COLUMNS : Slit events so the amr annotation event can be found later
        event_count = 1 
        for j in features_dict[ID]['annotations']:
            event_header = 'event_'+str(event_count)
            event_count += 1
            features_dict[ID][event_header] = [s for s in j if re.match(pattern_id, str(s))][0]
    count+=1
    
features_df = pd.DataFrame.from_dict(features_dict,orient='index')
features_df = features_df.drop(['location','protein_translation'],axis=1)
features_df = features_df.rename(columns={"id": "feature_id"})

# CREATE ANNOTATION EVENT DATAFRAME - LIST OF ANNOTAITON TOOLS (AKA EVENTS), EXPLUDE NON IMPORTANT TOOLS LIKE PRODIGAL
analysis_events_dict = {}
for i in annotation_genome['analysis_events']:
    analysis_events_dict[i['id']]= i 
analysis_events_df = pd.DataFrame.from_dict(analysis_events_dict, orient='index')

# AS THE EVENTS WERE SPLIT (LINE 170) INTO DIFFERENT COLUMNS, WE NEED TO DECIDE WHICH EVENT DID THE AMR ANNOTATION 
# collect ALL of the events that were involved with the AMR hits, 
collect_events = []

for i in features_df.columns:
    if "event_" in i:
        collect_events = collect_events + features_df[i].unique().tolist()
# make a unique list        
collect_events = list(set(collect_events))
# remove nan
collect_events = [x for x in collect_events if x == x]

# Make a subset ataframe of the events that are associated with AMR hits. 
for i in collect_events:
    analysis_events = analysis_events_df.loc[collect_events]

# Remove these rows that have unnessesary tools, we are left with list of tools that annotated the AMR
remove_these_events = ['GenomeAnnotation::renumber_features', 'glimmer3','prodigal']
# other_analysis_events = analysis_events[analysis_events['tool_name'].isin(remove_these_events)].unique().tolist()
other_analysis_events = analysis_events[analysis_events['tool_name'].isin(remove_these_events)]
# Collect all the column names that start with 'event_', there can be multiple. We need this information because we will merge them and decide what event annotated the AMR
event_collection = []
for i in features_df.columns: 
    if "event_" in i:
        event_collection.append(i)

#Replace all of the events that are in other_analysis_events and in the 'event_' columns wiht nan
features_df[event_collection] = features_df[event_collection].replace(other_analysis_events['id'].unique().tolist(), np.nan)

# Function to merge 'event_' columns based on conditions
def merge_columns(row, events):
    if events ==3:
        # Filter out the non-NaN values
        non_empty_values = [val for val in [row['event_1'], row['event_2'], row['event_3']] if pd.notna(val)]
        
        if len(non_empty_values) == 0:
            return np.nan  # All values are NaN
        elif len(set(non_empty_values)) == 1:
            return non_empty_values[0]  # All values are the same
        else:
            return np.random.choice(non_empty_values)  # Randomly pick a value if different
    if events ==4:
        # Filter out the non-NaN values
        non_empty_values = [val for val in [row['event_1'], row['event_2'], row['event_3'], row['event_4']] if pd.notna(val)]
        
        if len(non_empty_values) == 0:
            return np.nan  # All values are NaN
        elif len(set(non_empty_values)) == 1:
            return non_empty_values[0]  # All values are the same
        else:
            return np.random.choice(non_empty_values)  # Randomly pick a value if different

# Apply the function to each row in the DataFrame
features_df['merged_event'] = features_df.apply(merge_columns, args=(len(event_collection),) ,axis=1)
features_df = features_df.drop(columns=event_collection)


# MAP THE AMR TO THE ANNOTATION FEATURE EVENT AND TO THE EVENT TOOL
mapped = quality_json_df.merge(right=features_df, right_on='feature_id', left_on='feature_id',how='left')
mapped = mapped.merge(right=analysis_events_df, right_on='id', left_on='merged_event',how='left')
mapped.to_csv('./bvbrc_amr_annotation.tsv', sep='\t')

# MAKE CSV OF ADABOOST PREDICTED RESISTANCE
with open(amr_json, 'r') as file:
    resistance = json.load(file)
resistance_df = pd.DataFrame.from_dict(resistance)
resistance_df.to_csv('./bvbrc_predicted_resistance.tsv',sep='\t')