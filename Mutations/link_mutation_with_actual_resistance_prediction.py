#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np


# In[2]:


location_in_json = pd.read_csv('/home/kin672/gentb-summer22/Mutations/GenTB_Important_Mutation_Locations.csv')
impt_mutations = pd.read_csv('/home/kin672/gentb-summer22/Mutations/GenTB_Important_Mutations_Only.csv')
mutation_to_WHO = pd.read_csv('/home/kin672/gentb-summer22/Mutations/GenTB_Important_Mutations_to_WHOList.csv')
resistance = pd.read_csv('/home/kin672/gentb-summer22/Creating Summaries from Relevant Files/7-18 Output/full_df_7-18.csv')


# In[4]:


# End goal columns: Isolate, Mutation, Drug, Confidence, Resistant to Drug?, Index of Sublist, Index within Sublist
target = mutation_to_WHO[['Isolate', 'GenTB_Mutation', 'Drug', 'Confidence']]


# In[3]:


def resistant_to_drug(isolate, drug, df = resistance):
    resistance = pd.read_csv('/home/kin672/gentb-summer22/Creating Summaries from Relevant Files/7-18 Output/full_df_7-18.csv')
    if pd.isna(drug):
        return('N/A')
    try:
        return(df.loc[(df.ID == isolate) & (df.Drug == drug), 'Resistant'].values[0])
    except:
        print('Isolate: ' + str(isolate) + ', Drug: ' + str(drug))
        return('Problem')
    
# Exclude if no lineage:
exclude = ['ERR2677407', 'por1A', 'X101_S12', 'M45_S6', 'K67_S2', 'ERR2677386', 'kr-s1']

resistance_col = []
for index, row in target.iterrows():
    isolate = row['Isolate']
    if isolate in exclude:
        continue
    drug = row['Drug']
    resistance_col.append(resistant_to_drug(isolate, drug))
target['Resistant_to_Drug'] = resistance_col

# In[11]:


def sublist_indices(isolate, mutation, df = location_in_json):
    try:
        sublist_index = df.loc[(df.Isolate == isolate) & (df.Mutation == mutation), 'Sublist_index'].values[0]
        index_in_sublist = df.loc[(df.Isolate == isolate) & (df.Mutation == mutation), 'Index_within_sublist'].values[0]
    except:
        return([mutation, mutation])
    return([sublist_index, index_in_sublist])

sublist_col = []
index_in_sublist_col = []
for index, row in target.iterrows():
    isolate = row['Isolate']
    mutation = row['GenTB_Mutation']
    indices = sublist_indices(isolate, mutation)
    sublist_col.append(indices[0])
    index_in_sublist_col.append(indices[1])
    
target['Index_of_Sublist'] = sublist_col
target['Index_within_Sublist'] = index_in_sublist_col


# In[ ]:


target.to_csv('/home/kin672/gentb-summer22/Mutations/GenTB_Important_Mutation_with_Predicted_Resistance.csv')

