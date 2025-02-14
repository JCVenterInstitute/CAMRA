import pandas as pd
import logging
from itertools import combinations

logger = logging.getLogger("term_consolidation")

# Step 4 - Generate metadata on tools and databases
def generate_metadata(hamr_output_df):
    """Generate the tool-db metadata and display results."""
    frequency_df = find_frequencies(hamr_output_df)
    matrix_df = create_matching_hit_matix(hamr_output_df)
    logger.info(f"> Number of Hits found by each software-database permutation: \n{frequency_df.to_markdown()}")
    logger.info(f"> Concurring Software-Database Matrix.\n{matrix_df.to_markdown()}")
    return frequency_df, matrix_df

def create_matching_hit_matix (df: pd.DataFrame) -> pd.DataFrame:
    '''
    Creates a matrix that counts how often each pair of tool-database combinations 
    is observed within the same term consolidation group. The result can be used 
    for visualization (e.g., heatmaps).
    '''
    
    # For each loci group all the software/db associated to the loci. 
    grouped : pd.Series = df.groupby('loci_groups')['permutations'].apply(list)
    
    # Generate all combinations of tool-database pairs for each group
    pair_counts : dict = {}
    for pairs in grouped:
        unique_pairs = list(set(pairs))  # Remove duplicates in the same group
        for pair1, pair2 in combinations(unique_pairs,2):  # Find all the combinations of unique_pairs
            key : tuple = tuple(sorted([pair1, pair2]))  # Ensure order-independent keys # Make each combination a key 
            pair_counts[key] : int = pair_counts.get(key, 0) + 1 # For this key in the pair_count dictionary add +1 every time a loci exhibits this combination
    
    # Convert pair_counts dictionary to a DataFrame
    pair_counts_df : pd.DataFrame = pd.DataFrame(pair_counts.items(), columns=['pair', 'count'])

    # Make Heatmap
    # Split pairs into separate columns for clarity (optional)
    pair_counts_df[['tool_db_1', 'tool_db_2']] = pd.DataFrame(pair_counts_df['pair'].tolist(), index=pair_counts_df.index)
    pair_counts_df.drop(columns='pair', inplace=True)
    
    # Reshape into a matrix format (pivot)
    matrix_df = pair_counts_df.pivot(index='tool_db_1', columns='tool_db_2', values='count').fillna(0)
    return matrix_df  

def find_frequencies (df: pd.DataFrame) -> pd.DataFrame:   
    
    '''
    Analyzes the frequency of each tool-database combination (permutation) and 
    the frequency of unique hits within these combinations.
    '''
    
    def count_hit_frequency_of_software_db_permutations (df: pd.DataFrame) -> pd.DataFrame:
        # Create new 'frequency' df that value counts the permutations of tuple(software,database). 
        frequency_df : pd.DataFrame = df['permutations'].value_counts().reset_index()
        frequency_df.columns = ['permutation', 'permutation_frequency']
        return frequency_df
        
    def count_unique_hit_frequency_of_software_db_permutations (df: pd.DataFrame) -> pd.DataFrame:
        # Identify the frequency a loci group was found by a tuple(software,database)
        group_sizes : pd.Series = df.groupby('loci_groups').size()
        # Store the the groups that have a size == 1.
        onehit_loci : pd.Index = group_sizes[group_sizes == 1].index 
        # Create new 'single_groups' df that contains the rows that belong to the onehit_loci
        onehit_loci_df = df[df['loci_groups'].isin(onehit_loci)].copy()

        # Create new 'frequency' df that value counts the permutations of tuple(software,database). 
        frequency_df = onehit_loci_df['permutations'].value_counts().reset_index()
        frequency_df.columns = ['permutation', 'unique_hit_frequency']
        return frequency_df
    
    def merge_frequency_df (all_hits_df: pd.DataFrame, single_hit_df: pd.DataFrame) -> pd.DataFrame:
        # Merge the two dataframes on the 'permutation' column
        combined_df = pd.merge(
            all_hits_df,
            single_hit_df,
            on='permutation',
            how='left'
        ).fillna(0)  # Fill NaN for unique hits with 0 (no unique hits for some permutations)
                # Ensure frequencies are integers
        combined_df['unique_hit_frequency'] = combined_df['unique_hit_frequency'].astype(int)
        combined_df['permutation_frequency'] = combined_df['permutation_frequency'].astype(int)
        return combined_df
        
    def create_column (df: pd.DataFrame, new_column: str, func):
        # Create a new column 'permunations'
        df[new_column] : pd.Series = df.apply(
            func,
            axis=1)
        return df 
    
    if 'permutations' not in df.columns:
        create_column(df, 'permutations', lambda row: tuple(
                sorted([    
                    row['analysis_software_name'], 
                    row['reference_database_name']
                ])))

    software_db_frequency : pd.DataFrame = count_hit_frequency_of_software_db_permutations(df)    
    software_db_frequency_uniquehit : pd.DataFrame = count_unique_hit_frequency_of_software_db_permutations (df)    
    frequency_df : pd.DataFrame = merge_frequency_df(software_db_frequency, software_db_frequency_uniquehit)
    return frequency_df

def match_metadata (df:pd.DataFrame, column:str):
    return df[column].value_counts(dropna=False).to_markdown()