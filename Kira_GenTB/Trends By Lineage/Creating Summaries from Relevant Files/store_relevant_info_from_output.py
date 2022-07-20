#!/usr/bin/env python
# coding: utf-8

# In[95]:


import pandas as pd
import numpy as np
import json
from pandas.api.types import CategoricalDtype


# In[2]:


with open('//home//kin672//Kira_GenTB//Isolating Relevant Files//relevant_file_names.json') as f:
    relevant_files = json.load(f) # This is a dictionary with format {STRAIN ID: [OUTPUT PATH, LINEAGE PATH]}

    


# In[50]:


prefix = '//n//groups//gentb_www//predictData//'

# RF probability thresholds (if probability is below threshold, it is susceptible to the drug)
thresholds = {'rif': 0.002, 'inh': 0.22, 'emb': 0.082,'pza': 0.013,'str':0.047, 'cap':0.25,'amk':0.6,'cip': 0.42,'kan': 0.63,'levo':0.41,'oflx':0.33,'pas': 0.001,'eth': 0.32}


# In[226]:
catch_problems2 = []
catch_problems = []
catch_problems3 = []
master_df = pd.DataFrame({'ID':[], 'Drug':[], 'Resistant':[], 'Lineage':[]})

for strainID in relevant_files:
    value = relevant_files.get(strainID)
    json_path = prefix + value['Folder'] + '//' + value['Resistance']
    lineage_path = prefix + value['Folder'] + '//' + value['Lineage']
        
    # Parse lineage
    try:
        lineage = pd.read_csv(lineage_path, sep = '\t')
    except:
        catch_problems3.append(lineage_path)
        continue
    if any(['freschi' in string for string in lineage.columns]):
        freschi_lineage = str(lineage.loc[0, list(lineage.columns)[list(np.where(['freschi' in string for string in lineage.columns])[0])[0]]])
        if '(1/1)' in freschi_lineage:
            freschi_lineage = freschi_lineage.replace('(1/1)','') 
    else:
        catch_problems3.append(lineage_path)
        freschi_lineage = 'unfilled'
    
    # Parse json
    with open(json_path) as f:
        resistance = json.load(f)
        
    # Cut out the extra output at the bottom
    resistance = resistance[0]
    if len(resistance) != 13:
        catch_problems2.append(json_path)
        continue
    
    # Pull binary resistance outcome and probability per drug!
    for drug_index in range(len(resistance)):
        profile = resistance[drug_index]
        drug_name = profile[1]
        resistant = 'unfilled'
        
        # If the json file has only three outputs per drug, we need to manually determine the binary output using thresholds.
        if len(profile) == 5: # 3 outputs - probability is at index 2
            if float(profile[2]) < thresholds[drug_name]:
                resistant = '0'
            else:
                resistant = '1'
        elif len(profile) == 6: # 4 outputs - binary output is at index 2
            resistant = profile[2] 
        else:
            catch_problems.append(strainID + '/' + drug_name)
            continue
        
        # Add a line to the dataframe
        master_df = pd.concat([master_df, pd.DataFrame({'ID':[strainID], 'Drug':[drug_name], 'Resistant':[resistant], 'Lineage':[freschi_lineage]})], ignore_index = True)
    
        




# In[230]:


print('Total number of strains processed: ' + str(len(master_df)) + '\n')
print('Doesnt have 3 or 4 outputs per drug: ' + ', '.join(catch_problems))
print('Too many drugs in json: ' + ', '.join(catch_problems2))
print('Messed up lineage files: ' + ', '.join(catch_problems3))


# In[232]:


drug_order = CategoricalDtype(thresholds.keys(), ordered = True)


# In[234]:


master_df['Drug'] = master_df['Drug'].astype(drug_order)
master_df = master_df.sort_values(by = ['Drug', 'Lineage','Resistant'])

master_df.to_csv('//home//kin672//Kira_GenTB//full_df_attempt3.csv')



# In[242]:


summarized = master_df.groupby(['Drug','Lineage','Resistant']).count()
summarized = summarized.rename(columns={'ID':'NumStrains'})
summarized.reset_index(inplace=True)


# In[258]:


reformat_summarized = pd.DataFrame({'Drug':[], 'Lineage':[], 'Percentage of Strains that are Resistant':[], 'Number of Resistant Strains':[], 'Number of Susceptible Strains':[]})
for drug in np.unique(summarized.Drug):
    for lineage in np.unique(summarized.Lineage):
        num_resistant = int(summarized[(summarized.Drug == drug) & (summarized.Lineage == lineage) & (summarized.Resistant == '1')].NumStrains)
        num_susceptible = int(summarized[(summarized.Drug == drug) & (summarized.Lineage == lineage) & (summarized.Resistant == '0')].NumStrains)
        temp = pd.DataFrame({'Drug':[drug], 'Lineage':[lineage], 'Percentage of Strains that are Resistant':[num_resistant * 100 / (num_resistant + num_susceptible)], 'Number of Resistant Strains':[str(num_resistant)], 'Number of Susceptible Strains':[str(num_susceptible)]})
        reformat_summarized = pd.concat([reformat_summarized, temp], ignore_index=True)
reformat_summarized


# In[ ]:


reformat_summarized.to_csv('//home//kin672//Kira_GenTB//summarized_full_attempt3.csv')

