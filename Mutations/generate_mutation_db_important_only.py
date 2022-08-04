# Extracting only the mutations deemed by GenTB to be important (The first dictionary in the JSON)

import pandas as pd
import numpy as np
import json
from pandas.api.types import CategoricalDtype
import re
from utils import *


#
# Load WHO Database
#
mutations = pd.read_csv('/n/data1/hms/dbmi/farhat/Sanjana/MIC_data/WHO_resistance_variants_all.csv')

# Switching the abbreviations to the ones that I've been using, and removing any that I haven't been working with
drug_translation = {'AMI': 'amk', 'BDQ': 'bdq', 'CAP':'cap', 'CFZ':'cfz', 'DLM':'dlm', 'EMB':'emb', 'ETH':'eth', 'INH':'inh', 'KAN':'kan',
       'LZD':'lzd', 'PZA':'pza', 'RIF':'rif', 'STM':'str', 'LEV':'levo', 'MXF': 'moxi'}
mutations.drug = [drug_translation[original] for original in mutations.drug]
mutations = mutations.loc[[drug in ['amk', 'cap', 'emb', 'eth', 'inh', 'kan', 'rif', 'pza', 'str', 'levo'] for drug in mutations.drug],:]
mutations['searchable_variant'] = [str(i).replace('-','').upper() for i in list(mutations.variant)]

#
# Load all JSON files
#
with open('//home//kin672//gentb-summer22//Isolating Relevant Files//7-18 Output//relevant_file_names_7-18.json') as f:
    relevant_files = json.load(f) # This is a dictionary with format {STRAIN ID: {Folder: , Resistance: , Lineage: }

problems = {'The JSON is weird': [], 'This mutation is weird': {'Isolate':[], 'Mutation':[]}}

prefix = '//n//groups//gentb_www//predictData//'

# Loop through every unique isolate to isolate the mutation list
concatenate = []
            
for strainID in relevant_files:
    value = relevant_files.get(strainID)
    json_path = prefix + value.get('Folder') + '//' + value.get('Resistance')
    
    # Parse json
    with open(json_path) as f:
        resistance = json.load(f) 
    
    # Get full list of IMPORTANT resistance-related mutations for this isolate
    try:
        mut_list = [i for sublist in resistance[1].values() for sublist2 in sublist for i in sublist2 if i != None and i != 'Null']
        mut_list = list(pd.Series(mut_list).unique())
    except:
        problems['The JSON is weird'].append(str(strainID))
        continue
    
    # Parse through the mutations in the list
    # Don't even attempt to parse if the list is empty
    if len(mut_list) == 0:
        continue
    
    # Use Martin's method to break down the mutation, put output in a dataframe
    for mut in mut_list:
        try:  
            dict_breakdown = match_snp_name(mut)
        except:
            problems['This mutation is weird']['Isolate'].append(str(strainID))
            problems['This mutation is weird']['Mutation'].append(str(mut))
            continue
        
        temp = pd.DataFrame([dict_breakdown])
        temp['Isolate'] = strainID
        temp['GenTB_Mutation'] = mut
            
        concatenate.append(temp)
        

# Tie it all together
final = pd.concat(concatenate)
            
            
     
connection_to_wholist = {'Isolate':[], 'GenTB_Mutation':[], 'Drug':[], 'Confidence':[], 'Gene':[], 'WHO_Variant':[]}            

for index, temp in final.iterrows():
    candidates = mutations.loc[mutations.genome_index == temp['ntpos'],:]

    # Make sure the gene aligns
    temp_gene = temp['gene']
    if pd.isna(temp_gene):
        temp_gene = temp['rgene']

    candidates = candidates.loc[[gene in temp_gene for gene in candidates.gene], :]

    # Isolate the candidates who have one of the swaps indicated
    swaps = [temp['codes'], temp['coding'], temp['amino']]
    candidates = candidates.loc[[any(str(code) in str(searchable_variant) for code in swaps) for searchable_variant in candidates.searchable_variant], :]  

    # Move on if there are no relevant mutations
    if len(candidates) == 0:
        connection_to_wholist['Isolate'].append(temp['Isolate'])
        connection_to_wholist['GenTB_Mutation'].append(temp['GenTB_Mutation'])
        connection_to_wholist['Drug'].append('N/A')
        connection_to_wholist['Confidence'].append('N/A')
        connection_to_wholist['Gene'].append('N/A')
        connection_to_wholist['WHO_Variant'].append('N/A')
    
    # Add all valid candidates to the list
    for index, candidate in candidates.iterrows():
        connection_to_wholist['Isolate'].append(temp['Isolate'])
        connection_to_wholist['GenTB_Mutation'].append(temp['GenTB_Mutation'])
        connection_to_wholist['Drug'].append(candidate['drug'])
        connection_to_wholist['Gene'].append(candidate['gene'])
        connection_to_wholist['Confidence'].append(candidate['confidence'])
        connection_to_wholist['WHO_Variant'].append(candidate['variant'])


## INCLUDES ALL MUTATIONS WITHOUT REGARD TO WHO LIST
final.to_csv('/home/kin672/gentb-summer22/Mutations/GenTB_Important_Mutations_Only.csv', index = False)
pd.DataFrame(problems['This mutation is weird']).to_csv('/home/kin672/gentb-summer22/Mutations/GenTB_Important_Weird_Mutations.csv', index = False)
pd.DataFrame(connection_to_wholist).to_csv('/home/kin672/gentb-summer22/Mutations/GenTB_Important_Mutations_to_WHOList.csv', index = False)
    