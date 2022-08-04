## Find the locations of the mutations so they can potentially be related back to which drug GenTB is intending them for

import pandas as pd
import numpy as np
import json
from pandas.api.types import CategoricalDtype
import re
from utils import *

#
# Load all JSON files
#
with open('//home//kin672//gentb-summer22//Isolating Relevant Files//7-18 Output//relevant_file_names_7-18.json') as f:
    relevant_files = json.load(f) # This is a dictionary with format {STRAIN ID: {Folder: , Resistance: , Lineage: }

prefix = '//n//groups//gentb_www//predictData//'


# Loop through every unique isolate
problems = []
mut_list = {'Isolate':[], 'Mutation':[], 'Sublist_index':[], 'Index_within_sublist':[]}   
for strainID in relevant_files:
    value = relevant_files.get(strainID)
    json_path = prefix + value.get('Folder') + '//' + value.get('Resistance')
    
    # Parse json
    with open(json_path) as f:
        resistance = json.load(f) 
     
    try:
        list_of_sublists = list(resistance[1].values())[0] # Length is always 1
        for sublist_index in range(len(list_of_sublists)):
            sublist_mut_list = list(list_of_sublists[sublist_index])
            for index_within_sublist in range(len(sublist_mut_list)):
                mutation = sublist_mut_list[index_within_sublist]
                if (str(mutation) not in ['None', 'nan']) & (not pd.isna(mutation)):
                    mut_list['Isolate'].append(strainID)
                    mut_list['Mutation'].append(mutation)
                    mut_list['Sublist_index'].append(sublist_index)
                    mut_list['Index_within_sublist'].append(index_within_sublist)
    except:
        problems.append(str(strainID))
        continue
        
print(problems)
    
pd.DataFrame(mut_list).to_csv('/home/kin672/gentb-summer22/Mutations/GenTB_Important_Mutation_Locations.csv', index = False)
            