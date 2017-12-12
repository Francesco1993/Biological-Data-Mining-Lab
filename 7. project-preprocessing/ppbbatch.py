"""
PreProcessingBiological Batch Module
"""

import GEOparse
import pandas as pd
import numpy as np
import os
import unicodedata

# HELPERS
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
 
    return False

def load_dataset(dataset_id, download_location="."):
    """
    Load the dataset from disk (or download it if it does not exists)
    Arguments:
    - dataset_id: the ID of the dataset to load
    
    Output:
    - GSE object (GEOparse Series)
    """
    path = download_location + "/" + dataset_id + "_family.soft.gz"
    if os.path.exists(path):
        # Load from an existing file
        print("- Loading from", path)
        gse = GEOparse.get_GEO(filepath=path)
    else:
        # Download GSE and load it
        print("- Downloading", dataset_id)
        gse = GEOparse.get_GEO(geo=dataset_id, destdir=download_location + "/")
    return gse

def extract_dataframes(gse):
    """
    Returns a list of pandas dataframes
    once for "sample" (person) in the dataset
    """
    return [gsm.table for gsm_name, gsm in gse.gsms.items()]

def create_david_input(dataframes, column_name, destlocation='./input_for_david.txt'):
    """
    Create a txt of unique probes_id for david (online web app converter to entrez_id)
    """
    list_of_probes = []
    for dataframe in dataframes:
        list_of_probes += list(dataframe[column_name].map(str))

    set_of_unique_probes = sorted(set(list_of_probes))

    with open(destlocation, 'w') as f:
        for probe_id in set_of_unique_probes:
            f.write(probe_id + '\n')

    print('saved on', destlocation)

def create_mapper_from_david(david_output_path):
    """
    Returns a python dictionary that represents our mapper object
    Important: not all probs_id are mapped to an ENTREZ_GENE_ID
    probs_id without an enterez_id are not added to the dictionary
    """

    david_conversion_dataframe = pd.read_csv(david_output_path, sep = "\t")
    david_conversion_dataframe.columns = ["ID", "ENTREZ_GENE_ID", "Species", "Gene Name"]

    mapper = {}
    for index, row in david_conversion_dataframe.iterrows():
        probs_id = str(row['ID'])

        if not row[destination_label] or pd.isnull(row[destination_label]):
            continue
        
        # then convert to string
        if not is_number(str(row['ENTREZ_GENE_ID'])):
            #print('error', str(row[destination_label]))
            continue
        else:
            new_value = str(int(row['ENTREZ_GENE_ID']))

        if probs_id in mapper and mapper[probs_id] != new_value:
            # Multiple enterez id for the same probs
            # Set their value to None to invalid them
            # Elements set to "None" are then removed
            mapper[probs_id] = None
            
        if probs_id not in mapper and not pd.isnull(row[destination_label]):
            mapper[probs_id] = new_value

    # Remove invalid mapping (value = None)
    # (Some of the probes are linked with multiple numbers (enterez_id ?) using /// as separator)
    filtered_mapper = {k:v for k,v in mapper.items() if v != None and '/' not in v}
    
    return filtered_mapper

def create_mapper_from_platform(gse, from_label, destination_label):
    """   
    Returns a python dictionary that represents our mapper object
    Important: not all probs_id are mapped to an ENTREZ_GENE_ID
    probs_id without an enterez_id are not added to the dictionary
    """

    meta_data_tables = []
    for gpl_name, gpl in gse.gpls.items():
        meta_data_tables.append(gpl.table)

    mapper = {}
    for df in meta_data_tables:
        for index, row in df.iterrows():
            # convert probes_id to string
            probs_id = str(row[from_label])

            if not row[destination_label] or pd.isnull(row[destination_label]):
                continue
            
            # then convert to string
            if not is_number(str(row[destination_label])):
                #print('error', str(row[destination_label]))
                continue
            else:
                new_value = str(int(row[destination_label]))
            
            if probs_id in mapper and mapper[probs_id] != new_value:
                # Multiple enterez id for the same probs
                # Set their value to None to invalid them
                # Elements set to "None" are then removed
                mapper[probs_id] = None
                
            if probs_id not in mapper and not pd.isnull(row[destination_label]):
                mapper[probs_id] = new_value
            
    # Remove invalid mapping (value = None)
    # (Some of the probes are linked with multiple numbers (enterez_id ?) using /// as separator)
    filtered_mapper = {k:v for k,v in mapper.items() if v != None and '/' not in v}
    
    return filtered_mapper

def mapper_to_pandas_df(mapper, from_label, destination_label):
    """
    Convert mapper to pandas dataframe
    """
    mapper_df = pd.DataFrame.from_dict(mapper, orient='index')
    mapper_df.index.name = from_label
    mapper_df.columns = [destination_label]  
    return mapper_df

def float_or_nan(x):
    if x == np.nan or str(x).lower() == 'null':
        return np.nan
    else:
        return float(x)

def remove_inf(x):
    if np.isinf(x):
        return np.nan
    else:
        return float(x)

def apply_log2(x):
    return np.log2(x)

def filter_and_normalize(data_frames, mapper, probs_label, sample_indexes, isLog2=True):
    data_frames_entrez = []
    for old_df in data_frames:
        df = old_df.filter(items=[probs_label, 'VALUE'])
        df['VALUE'] = df['VALUE'].apply(float_or_nan)
        df['VALUE'] = df['VALUE'].apply(remove_inf)
        df = df.dropna(axis=0, how='any')
        df[probs_label] = df[probs_label].map(str)
        df = pd.merge(df, mapper, how='inner', left_on=[probs_label], right_index=True, sort=False)

        if isLog2:
            df['VALUE'] = 2**df['VALUE']

        df = df.groupby('ENTREZ_GENE_ID').mean()

        if np.isinf(df.VALUE.sum()):
            print(df.VALUE)
        
        df.VALUE = df.VALUE / df.VALUE.sum()
        
        df['VALUE'] = df['VALUE'].apply(apply_log2)
        
        df['VALUE'] = df['VALUE'].apply(remove_inf)
        df = df.dropna(axis=0, how='any')
    
        data_frames_entrez.append(df)

    # concat data_frames by columns and use the people's label as index
    merged_entrez_value_df = pd.concat(data_frames_entrez, axis=1, keys=sample_indexes)

    # remove gene with at least 1 missing value
    merged_entrez_value_df = merged_entrez_value_df.dropna(axis=0, how='any')
    merged_entrez_value_df.columns = merged_entrez_value_df.columns.droplevel([1])
    return merged_entrez_value_df