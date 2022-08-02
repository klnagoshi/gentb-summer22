#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd
import numpy as np


# In[47]:


all_mutations = pd.read_csv('/home/kin672/gentb-summer22/Mutations/all_mutations.csv')
mutation_count = pd.read_csv('/home/kin672/gentb-summer22/Mutations/Mutation_Count.csv')
full_db = pd.read_csv('/home/kin672/gentb-summer22/Creating Summaries from Relevant Files/7-18 Output/full_df_7-18.csv')


# In[ ]:


# Summarize mutation identities for each isolate
num_mutations_per_drug = {'Isolate':[], 'total':[], 'rif':[], 'inh':[], 'emb':[], 'pza':[], 'str':[], 'cap':[], 'amk':[], 'cip':[], 'kan':[], 'levo':[], 'oflx':[], 'pas':[], 'eth':[], 'unknown':[]}

for isolate in full_db.ID.unique():
    num_mutations_per_drug['Isolate'].append(isolate)
    if isolate in mutation_count.Isolate.unique():
        num_total_mutations = mutation_count.loc[mutation_count['Isolate'] == isolate, 'Num_Mutations'].values[0]
        rif = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'rif')])
        inh = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'inh')])
        emb = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'emb')])
        pza = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'pza')])
        strep = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'str')])
        cap = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'cap')])
        amk = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'amk')])
        cip = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'cip')])
        kan = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'kan')])
        levo = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'levo')])
        oflx = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'oflx')])
        pas = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'pas')])
        eth = len(all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == 'eth')])
        not_specified = num_total_mutations - len(all_mutations.loc[all_mutations['isolate'] == isolate, 'GenTB Mutation'].values.unique()) 
        total = num_total_mutations
    else:
        total, rif, inh, emb, pza, strep, cap, amk, cip, kan, levo, oflx, pas, eth, not_specified = 0
    
    num_mutations_per_drug['total'].append(total)
    num_mutations_per_drug['rif'].append(rif)
    num_mutations_per_drug['inh'].append(inh)
    num_mutations_per_drug['emb'].append(emb)
    num_mutations_per_drug['pza'].append(pza)
    num_mutations_per_drug['str'].append(strep)
    num_mutations_per_drug['cap'].append(cap)
    num_mutations_per_drug['amk'].append(amk)
    num_mutations_per_drug['cip'].append(cip)
    num_mutations_per_drug['kan'].append(kan)
    num_mutations_per_drug['levo'].append(levo)
    num_mutations_per_drug['oflx'].append(oflx)
    num_mutations_per_drug['pas'].append(pas)
    num_mutations_per_drug['eth'].append(eth)
    num_mutations_per_drug['unknown'].append(not_specified)
    
pd.DataFrame(num_mutations_per_drug).to_csv('/home/kin672/gentb-summer22/Mutations/Num_Mutations_Per_Drug.csv', index = False)

    


# In[51]:


# Function that produces confidence db for a single drug
# Product db includes (for each isolate): Resistance, Lineage, # of WHO-identified mutations, # of high confidence mutations, # low confidence
def summarize_mutations(drug, all_mutations = all_mutations, mutation_count = mutation_count):
    db = {'Isolate':[], 'Lineage':[], 'Resistant':[], 'WHO-Identified':[], 'High_Confidence':[], 'Low_Confidence':[]}
    
    for isolate in full_db.ID.unique():
        lineage = full_db.loc[(full_db['ID'] == isolate) & (full_db['Drug'] == 'rif'), 'Lineage'].values[0]
        resistant = full_db.loc[(full_db['ID'] == isolate) & (full_db['Drug'] == 'rif'), 'Resistant'].values[0]
        db['Isolate'].append(isolate)
        db['Lineage'].append(lineage)
        db['Resistant'].append(resistant)
        
        # If Isolate has any WHO-identified mutations, classify them
        if isolate in all_mutations.isolate.unique():
            WHO_mutations = all_mutations[(all_mutations['isolate'] == isolate) & (all_mutations['drug'] == drug)]
            num_WHO_mutations = len(WHO_mutations)
            num_high_confidence = len(WHO_mutations[(WHO_mutations['confidence'] == '1) Assoc w R') | (WHO_mutations['confidence'] == '2) Assoc w R - Interim')])
            num_low_confidence = num_WHO_mutations - num_high_confidence
            
            db['WHO-Identified'].append(num_WHO_mutations)
            db['High_Confidence'].append(num_high_confidence)
            db['Low_Confidence'].append(num_low_confidence)
        else:
            num_WHO_mutations = 0
            num_high_confidence = 0 
            num_low_confidence = 0
            db['WHO-Identified'].append(num_WHO_mutations)
            db['High_Confidence'].append(num_high_confidence)
            db['Low_Confidence'].append(num_low_confidence)
            
    return(pd.DataFrame(db))
    


# In[ ]:


for drug in full_db.Drug.unique():
    pd.DataFrame(summarize_mutations(drug)).to_csv('/home/kin672/gentb-summer22/Mutations/Breakdown per Drug/' + str(drug) + '_summary_of_mutations.csv', index = False)

