#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import json
import os


# In[ ]:


# Create dictionary with format {STRAIN ID: [lineage.txt file, matrix.json file]}
target_dict = {}
IDs_with_errors = {'No lineage files':{}, 'More than one lineage file':{}, 'The strainID already existed':{}}
predictData = os.listdir('//n//groups//gentb_www//predictData')
empty_json = []

# Walk through every submission folder in predictData
for folder in predictData:
    files = os.listdir('//n//groups//gentb_www//predictData//' + folder)
    
    # Iterate through each file in each submission folder to pull relevant .json and .txt files only
    for file in files:
        if 'matrix.json' in file:
            json_path = '//n//groups//gentb_www//predictData//' + folder + '//' + file
            if open(json_path).read() == '':
                empty_json.append(folder + '/' + file)
                continue
            else:
                with open(json_path) as f:
                    resistance = json.load(f)
                strainID = resistance[0][0][0].split('/')[-1]
            
            # Find corresponding lineage.txt file 
            lineage_files = list(filter(lambda x: (strainID in x and 'lineage.txt' in x), files))
            
            # If no lineage file exists or too many exist, there is a problem.
            if len(lineage_files) == 0:
                if folder not in IDs_with_errors['No lineage files'].keys():
                    IDs_with_errors['No lineage files'][folder] = [strainID]
                elif strainID not in IDs_with_errors['No lineage files'].get(folder):
                    IDs_with_errors['No lineage files'][folder].append(strainID)
                continue
            elif len(lineage_files) != 1:
                if folder not in IDs_with_errors['More than one lineage file'].keys():
                    IDs_with_errors['More than one lineage file'][folder] = [strainID]
                elif strainID not in IDs_with_errors['More than one lineage file'].get(folder):
                    IDs_with_errors['More than one lineage file'][folder].append(strainID)
                continue
            
            # Add strainID to dictionary if not already in -- if it is already present, there is a problem.
            if strainID not in target_dict.keys():
                target_dict[strainID] = {'Folder':folder, 'Resistance':file, 'Lineage':lineage_files[0]}
            else:
                if folder not in IDs_with_errors['The strainID already existed'].keys():
                    IDs_with_errors['The strainID already existed'][folder] = [strainID]
                elif strainID not in IDs_with_errors['The strainID already existed'].get(folder):
                    IDs_with_errors['The strainID already existed'][folder].append(strainID)
            
            
IDs_with_errors['.json file is empty'] = empty_json


# In[ ]:


with open('//home//kin672//Kira_GenTB//Isolating Relevant Files//relevant_file_names.json', 'w') as f1:
    json.dump(target_dict, f1,  indent=4)
    
with open('//home//kin672//Kira_GenTB//Isolating Relevant Files//strainIDs_with_errors.json', 'w') as f2:
    json.dump(IDs_with_errors, f2,  indent=4)  

