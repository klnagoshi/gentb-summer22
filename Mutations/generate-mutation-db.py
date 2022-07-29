#!/usr/bin/env python
# coding: utf-8

# In[7]:


import pandas as pd
import numpy as np
import json
from pandas.api.types import CategoricalDtype


# In[8]:


mutations = pd.read_csv('/n/data1/hms/dbmi/farhat/Sanjana/MIC_data/WHO_resistance_variants_all.csv')


# In[9]:


# Switching the abbreviations to the ones that I've been using, and removing any that I haven't been working with
drug_translation = {'AMI': 'amk', 'BDQ': 'bdq', 'CAP':'cap', 'CFZ':'cfz', 'DLM':'dlm', 'EMB':'emb', 'ETH':'eth', 'INH':'inh', 'KAN':'kan',
       'LZD':'lzd', 'PZA':'pza', 'RIF':'rif', 'STM':'str', 'LEV':'levo', 'MXF': 'moxi'}
mutations.drug = [drug_translation[original] for original in mutations.drug]
mutations = mutations.loc[[drug in ['amk', 'cap', 'emb', 'eth', 'inh', 'kan', 'rif', 'pza', 'str', 'levo'] for drug in mutations.drug],:]
mutations['searchable_variant'] = [str(i).replace('-','').upper() for i in list(mutations.variant)]


# In[11]:


with open('//home//kin672//gentb-summer22//Isolating Relevant Files//7-18 Output//relevant_file_names_7-18.json') as f:
    relevant_files = json.load(f) # This is a dictionary with format {STRAIN ID: {Folder: , Resistance: , Lineage: }


# In[152]:


def break_down(string):
    broken = {'original': string, 'gene': [], '1': [], '2': [], 'genome_index': [], 'change': []}
    temp = string.split('_')
    broken['gene'] = temp.pop(-1)
    broken['1'] = temp.pop(0)
    broken['2'] = temp.pop(0)
    broken['genome_index'] = temp.pop(0)
    broken['change'] = temp

    return(broken)


# In[189]:


final_db = pd.DataFrame({'isolate':[], 'drug':[], 'confidence':[], 'gene':[], 'genome_index':[], 'GenTB Mutation':[], 'WHO Variant':[]})


# In[150]:


prefix = '//n//groups//gentb_www//predictData//'

# Loop through every unique isolate
for strainID in relevant_files:
    value = relevant_files.get(strainID)
    json_path = prefix + value.get('Folder') + '//' + value.get('Resistance')
    
    # Parse json
    with open(json_path) as f:
        resistance = json.load(f) 
    
    # Get full list of resistance-related mutations for this isolate
    mut_list = [i for sublist in resistance[1].values() for sublist2 in sublist for i in sublist2 if i != None]
    mut_list = mut_list + [i for sublist in resistance[2].values() for sublist2 in sublist for i in sublist2 if i != None]
    mut_list = mut_list + [item for sublist in resistance[3].values() for item in sublist if item != 'Null']
    mut_list = list(pd.Series(mut_list).unique())
    
    # Match the mutation with the corresponding mutation(s) in the WHO database
    for mut in mut_list:
        x = break_down(mut)
        candidates = mutations.loc[mutations.genome_index == x['genome_index'],:]
        
        # Make sure the gene aligns
        candidates = candidates.loc[[gene in x['gene'] for gene in candidates.gene], :]
        
        # Isolate the candidates who have one of the amino acid swaps indicated
        candidates = candidates.loc[[any(code in searchable_variant for code in x['change']) for searchable_variant in candidates.searchable_variant], :]  
        
        # Move on if there are no relevant mutations
        if len(candidates) == 0:
            continue
            
        temp = candidates.loc[:,['drug', 'genome_index', 'confidence', 'gene']]
        temp['WHO Variant'] = candidates.variant
        temp['isolate'] = strainID
        temp['GenTB Mutation'] = x['original']
        
        final_db = pd.concat(final_db, temp)
            
        
        
        
        


# In[190]:


final_db.to_csv('/home/kin672/gentb-summer22/Mutations/all_mutations.csv', index = False)

