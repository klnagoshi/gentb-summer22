#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import json
from pandas.api.types import CategoricalDtype

master_df = pd.read_csv('/home/kin672/gentb-summer22/Summarize Sanjana Dataset/full_control_db_7-25.csv')
master_df = master_df[[str(i)[0:1] in ['1','2','3','4','5','6','B'] for i in master_df.Lineage]]

master_df = master_df[[not any(pd.isna(master_df.loc[master_df.ID == i, 'Resistant'])) for i in master_df['ID']]]
print('Number of IDs included: ' + str(len(master_df['ID'].unique())))

summarized = master_df.groupby(['Drug','Lineage','Resistant']).count()
summarized = summarized.rename(columns={'ID':'NumStrains'})
summarized.reset_index(inplace=True)


# In[258]:


reformat_summarized = pd.DataFrame({'Drug':[], 'Lineage':[], 'Percentage of Strains that are Resistant':[], 'Number of Resistant Strains':[], 'Number of Susceptible Strains':[]})
for drug in np.unique(summarized.Drug):
    for lineage in np.unique(summarized.Lineage):
        num_resistant = list(summarized.loc[(summarized.Drug == drug) & (summarized.Lineage == lineage) & (summarized.Resistant == 'R'), 'NumStrains'])
        if (len(num_resistant) == 1):
            num_resistant = int(num_resistant[0])
        elif (len(num_resistant) == 0):
            num_resistant = 0
        else:
            num_resistant = 'Messed up'
        num_susceptible = summarized[(summarized.Drug == drug) & (summarized.Lineage == lineage) & (summarized.Resistant == 'S')].NumStrains
        if (len(num_susceptible) == 1):
            num_susceptible = int(num_susceptible[0])
        elif (len(num_susceptible) == 0):
            num_susceptible = 0
        else:
            num_susceptible = 'Messed up'
        temp = pd.DataFrame({'Drug':[drug], 'Lineage':[lineage], 'Percentage of Strains that are Resistant':[num_resistant * 100 / (num_resistant + num_susceptible)], 'Number of Resistant Strains':[str(num_resistant)], 'Number of Susceptible Strains':[str(num_susceptible)]})
        reformat_summarized = pd.concat([reformat_summarized, temp], ignore_index=True)


# In[ ]:


reformat_summarized.to_csv('/home/kin672/gentb-summer22/Summarize Sanjana Dataset/summarized_control_no_nas_8-1.csv', index = False)

