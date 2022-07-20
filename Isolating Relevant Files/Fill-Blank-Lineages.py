#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import json
import glob, os, subprocess


# In[36]:


without_blanks = pd.read_csv('/home/kin672/Kira_GenTB/Isolating Relevant Files/vcf_to_lineages.csv')
full_df = pd.read_csv('/home/kin672/Kira_GenTB/Creating Summaries from Relevant Files/full_df_7-14.csv')
with open('//home//kin672//Kira_GenTB//Isolating Relevant Files//relevant_file_names_v2.json') as f:
    file_dict = json.load(f)
    
file_dict


# In[33]:


blank_IDs = full_df.loc[pd.isna(full_df.Lineage), :].ID.unique()


# In[ ]:


for row in range(len(full_df)):
    if not pd.isna(full_df.at[row, 'Lineage']):
        continue
    else:
        ID = full_df.at[row, 'ID']
        temp = file_dict.get(ID)
        folder_path = '//n//groups//gentb_www//predictData//' + temp.get('Folder')
        files = os.listdir(folder_path)

        for file in files:
            freschi_lineage = 'BLANK'
            lineage_files = list(filter(lambda x: (ID in x and 'lineage.txt' in x), files))

            if len(lineage_files) == 1:
                lineage_path = folder_path + '//' + lineage_files[0]
                try:
                    lineage = pd.read_csv(lineage_path, sep = '\t')
                # If for some reason this lineage file is not readable, record it as a failure.
                except:
                    freschi_lineage = 'Cant do lineage caller or existing file'
                # Pull the freschi lineage if it exists
                if any(['freschi' in string for string in lineage.columns]):
                    freschi_lineage = str(lineage.loc[0, list(lineage.columns)[list(np.where(['freschi' in string for string in lineage.columns])[0])[0]]]).replace('(1/1)','') 
                else:
                    freschi_lineage = 'Cant do lineage caller or existing file'
            else:
                freschi_lineage = 'Cant do lineage caller or existing file'

            full_df.at[row, 'Lineage'] = freschi_lineage
            
        


# In[ ]:


pd.to_csv(full_df, '/home/kin672/Kira_GenTB/Creating Summaries from Relevant Files/full_df_7-18.csv', index = False)

