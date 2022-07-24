#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import json
from pandas.api.types import CategoricalDtype

master_df = pd.read_csv('/home/kin672/gentb-summer22/Summarize Sanjana Dataset/full_control_db_7-24.csv')

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


reformat_summarized.to_csv('/home/kin672/gentb-summer22/Summarize Sanjana Dataset/summarized_control_7-24.csv')

