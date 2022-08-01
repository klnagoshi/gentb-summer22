import pandas as pd
import numpy as np
import json
from pandas.api.types import CategoricalDtype
import re

mutations = pd.read_csv('/n/data1/hms/dbmi/farhat/Sanjana/MIC_data/WHO_resistance_variants_all.csv')

# Switching the abbreviations to the ones that I've been using, and removing any that I haven't been working with
drug_translation = {'AMI': 'amk', 'BDQ': 'bdq', 'CAP':'cap', 'CFZ':'cfz', 'DLM':'dlm', 'EMB':'emb', 'ETH':'eth', 'INH':'inh', 'KAN':'kan',
       'LZD':'lzd', 'PZA':'pza', 'RIF':'rif', 'STM':'str', 'LEV':'levo', 'MXF': 'moxi'}
mutations.drug = [drug_translation[original] for original in mutations.drug]
mutations = mutations.loc[[drug in ['amk', 'cap', 'emb', 'eth', 'inh', 'kan', 'rif', 'pza', 'str', 'levo'] for drug in mutations.drug],:]
mutations['searchable_variant'] = [str(i).replace('-','').upper() for i in list(mutations.variant)]

with open('//home//kin672//gentb-summer22//Isolating Relevant Files//7-18 Output//relevant_file_names_7-18.json') as f:
    relevant_files = json.load(f) # This is a dictionary with format {STRAIN ID: {Folder: , Resistance: , Lineage: }
    
# Function that breaks down a GenTB mutation into its component parts
def break_down(string):
    broken = {'original': string, 'gene': [], '1': [], '2': [], 'genome_index': [], 'change': []}
    temp = re.split('_|\\.', string)
    if not len(temp) >= 5:
        print(string)
    try:
        broken['gene'] = temp.pop(-1)
        broken['1'] = temp.pop(0)
        broken['2'] = temp.pop(0)
        broken['genome_index'] = temp.pop(0)
        broken['change'] = temp
    except:
        return(broken)

    return(broken)

# Set up dictionaries
final_db = pd.DataFrame({'isolate':[], 'drug':[], 'confidence':[], 'gene':[], 'genome_index':[], 'GenTB Mutation':[], 'WHO Variant':[]})
problems = {'The JSON is weird': [], 'This mutation is weird': {'Isolate':[], 'Mutation':[]}} # In theory, there should be no more weird JSONs
mutation_count = {'Isolate':[], 'Num_Mutations':[]}

prefix = '//n//groups//gentb_www//predictData//'

# Loop through every unique isolate
for strainID in relevant_files:
    value = relevant_files.get(strainID)
    json_path = prefix + value.get('Folder') + '//' + value.get('Resistance')
    
    # Parse json
    with open(json_path) as f:
        resistance = json.load(f) 
    
    # Get full list of resistance-related mutations for this isolate
    try:
        mut_list = [i for sublist in resistance[1].values() for sublist2 in sublist for i in sublist2 if i != None and i != 'Null']
        mut_list = mut_list + [i for sublist in resistance[2].values() for sublist2 in sublist for i in sublist2 if i != None and i != 'Null']
        mut_list = mut_list + [item for sublist in resistance[3].values() for item in sublist if item != None and item != 'Null']
        mut_list = list(pd.Series(mut_list).unique())
    except:
        try:
            mut_list = [i for sublist in resistance[1].values() for sublist2 in sublist for i in sublist2 if i != None and i != 'Null']
            mut_list = mut_list + [i for sublist in resistance[2].values() for sublist2 in sublist for i in sublist2 if i != None and i != 'Null']
        except:
            problems['The JSON is weird'].append(str(strainID))
            continue
    
    # Record the total number of mutations
    mutation_count['Isolate'].append(strainID)
    mutation_count['Num_Mutations'].append(len(mut_list))
    
    # Don't even attempt to parse if the list is empty
    if len(mut_list) == 0:
        continue
    
    # Match the mutation with the corresponding mutation(s) in the WHO database
    for mut in mut_list:
        x = break_down(mut)
        if x == 'Bad Breakdown' or [] in x.values():
            problems['This mutation is weird']['Isolate'].append(str(strainID))
            problems['This mutation is weird']['Mutation'].append(str(mut))
            continue
            
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
        
        final_db = pd.concat([final_db, temp], ignore_index = True)
        
        
final_db.to_csv('/home/kin672/gentb-summer22/Mutations/all_mutations.csv', index = False)
pd.DataFrame(problems['The JSON is weird']).to_csv('/home/kin672/gentb-summer22/Mutations/Weird_JSONS.csv', index = False)
pd.DataFrame(problems['This mutation is weird']).to_csv('/home/kin672/gentb-summer22/Mutations/Weird_Mutations.csv', index = False)
pd.DataFrame(mutation_count).to_csv('/home/kin672/gentb-summer22/Mutations/Mutation_Count.csv', index = False)
            