# Using Martin's GenTB default function to parse SNP names
# Checking if this is different from my str_split method

import pandas as pd
import numpy as np
import json
from pandas.api.types import CategoricalDtype
import re
from utils import *

with open('//home//kin672//gentb-summer22//Isolating Relevant Files//7-18 Output//relevant_file_names_7-18.json') as f:
    relevant_files = json.load(f) # This is a dictionary with format {STRAIN ID: {Folder: , Resistance: , Lineage: }

problems = {'The JSON is weird': [], 'This mutation is weird': {'Isolate':[], 'Mutation':[]}}

prefix = '//n//groups//gentb_www//predictData//'

# Loop through every unique isolate
concatenate = []

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
        temp['Full Mutation'] = mut
        concatenate.append(temp)

# Tie it all together
final = pd.concat(concatenate)




## INCLUDES ALL MUTATIONS WITHOUT REGARD TO WHO LIST
final.to_csv('/home/kin672/gentb-summer22/Mutations/Mutation_Info_Martin.csv', index = False)
pd.DataFrame(problems['This mutation is weird']).to_csv('/home/kin672/gentb-summer22/Mutations/Weird_Mutations_Martin.csv', index = False)
        
    