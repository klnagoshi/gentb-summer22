#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import json
import glob, os, subprocess


# In[6]:


# Read all files
target_cols = ['Isolate', 'RIFAMPICIN','ISONIAZID', 'ETHAMBUTOL', 'PYRAZINAMIDE', 'STREPTOMYCIN', 'CAPREOMYCIN', 'AMIKACIN', 'CIPROFLOXACIN', 'KANAMYCIN', 'LEVOFLOXACIN', 'OFLOXACIN', 'PARA-AMINOSALICYLIC_ACID', 'ETHIONAMIDE']
db = pd.read_csv('/n/data1/hms/dbmi/farhat/anna/focus_cnn/master_table_resistance.csv', usecols = target_cols)


# In[7]:


db = db.rename(columns={'Isolate':'Isolate', 'RIFAMPICIN':'rif','ISONIAZID':'inh', 'ETHAMBUTOL':'emb', 'PYRAZINAMIDE':'pza', 'STREPTOMYCIN':'str', 'CAPREOMYCIN':'cap', 'AMIKACIN':'amk', 'CIPROFLOXACIN':'cip', 'KANAMYCIN':'kan', 'LEVOFLOXACIN':'levo', 'OFLOXACIN':'oflx', 'PARA-AMINOSALICYLIC_ACID':'pas', 'ETHIONAMIDE':'eth'})
db.to_csv('/home/kin672/gentb-summer22/Summarize Sanjana Dataset/resistance_data.csv', index = False)


# In[48]:


# Fast Lineage Caller for all

## SANJANA'S LINEAGE CALLER METHOD, Edited
def get_lineages(fNames):
    
    isolates = []
    lineages = []
    
    for strainID in fNames:
        try:
            x = fNames.get(strainID)
            proc = subprocess.Popen(f"/home/kin672/anaconda3/envs/jupytervenv/bin/fast-lineage-caller {x} --noheader --count", shell=True, encoding='utf8', stdout=subprocess.PIPE)
            output = proc.communicate()[0]

            # the second value is the Freschi et al lineage
            freschi = output.split("\t")[1]
            freschi = freschi.replace('lineage', '')
            lineages.append(freschi)
            isolates.append(strainID)
        except:
            lineages.append('Lineage caller failed')
            isolates.append(strainID)
        
    return pd.DataFrame({"Isolate": isolates, "Lineage": lineages})


# In[83]:


y = pd.read_csv('/n/data1/hms/dbmi/farhat/rollingDB/genomic_data/strains_with_no_vcf', escapechar='\\', sep = '|', header = None)[0]
x = pd.read_csv('/n/data1/hms/dbmi/farhat/rollingDB/genomic_data/strains_with_no_vcf_no_errors', escapechar='\\', sep = '|', header = None)[0]
x = x[0:(len(x)-2)]
x


# In[93]:


no_vcfs = list(y) + list(x)
no_vcfs = [i.replace(' ', '') for i in no_vcfs]


# In[96]:


cryptic = os.listdir('/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/')
not_cryptic = os.listdir('/n/data1/hms/dbmi/farhat/rollingDB/genomic_data/')


# In[100]:


errors = {'Didnt have folder associated: ': [], 'Didnt have .vcf associated: ':[]}
vcfs = {}

for ID in db.Isolate:
    path = 'NONE'
    if ID not in no_vcfs:
        if ID in cryptic:
            path = '/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/' + str(ID) + '/pilon/' + str(ID) + '.vcf'
            vcfs[ID] = path
        elif ID in not_cryptic:
            path = '/n/data1/hms/dbmi/farhat/rollingDB/genomic_data/' + str(ID) + '/pilon/' + str(ID) + '.vcf'
            vcfs[ID] = path
        else: 
            errors['Didnt have folder associated: '].append(ID)
            continue
    else:
        errors['Didnt have .vcf associated: '].append(ID)
        continue
    
lineage_table = get_lineages(vcfs)


# In[8]:


master_df = pd.DataFrame({'ID':[], 'Drug':[], 'Resistant':[], 'Lineage':[]})

drugs = ['rif','inh','emb','pza','str','cap','amk','cip','kan','levo','oflx', 'pas','eth']

    


# In[ ]:


for ID in db.Isolate:
    
    # First, pull lineage value from the table generated earlier.
    if ID in errors.get('Didnt have folder associated: '):
        freschi_lineage = 'No folder'
    else:
        freschi_lineage = lineage_table.loc[lineage_table.Isolate == ID, 'Lineage'].to_list()
        freschi_lineage = [str(i).replace('(1/1)', '') for i in freschi_lineage]
        freschi_lineage = ', '.join(freschi_lineage)
    
    # Create a dataframe with all of the information for this isolate that will be added to the bottom of our master dataframe
    temp = pd.DataFrame(db.loc[db.Isolate == ID, drugs]).transpose().reset_index()
    temp.columns = ['Drug', 'Resistant']
    temp['Lineage'] = freschi_lineage
    temp.insert(0, 'ID', [ID]*13, allow_duplicates = True)
    
    master_df = pd.concat([master_df, temp])
    

# In[152]:


master_df.to_csv('/home/kin672/gentb-summer22/Summarize Sanjana Dataset/full_control_db_7-20.csv', index = False)

