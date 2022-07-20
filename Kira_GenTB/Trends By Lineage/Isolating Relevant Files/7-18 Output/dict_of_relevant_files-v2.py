#!/usr/bin/env python
# coding: utf-8

# In[106]:


#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import json
import glob, os, subprocess

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


# In[107]:



# Create dictionary with format {STRAIN ID: [lineage.txt file, matrix.json file]}
target_dict = {}
IDs_with_errors = {'No lineage files, no vcf':{}, 'More than one .vcf file, did not check for lineage files': {}, 'No vcf, More than one lineage file':{}, 'This strainID is present in more than one folder':{}}
predictData = os.listdir('//n//groups//gentb_www//predictData')
empty_json = []
vcf_dict = {}
not_rf = []

num_json = 0
unique_json = 0

test = []


# In[108]:

predictData.remove('access-filtered.log.gz')

# Walk through every submission folder in predictData
for folder in predictData:
    files = os.listdir('//n//groups//gentb_www//predictData//' + folder)
    
    # Iterate through each file in each submission folder to pull relevant .json and/or .txt files only
    for file in files:
        if 'matrix.json' in file:
            json_path = '//n//groups//gentb_www//predictData//' + folder + '//' + file
            vcf_check = False
            
            with open(json_path) as f:
                try:
                    resistance = json.load(f)
                except:
                    empty_json.append(folder + '/' + file)
                    continue
                    

            # Cut out the extra output at the bottom
            resistance = resistance[0]
            if len(resistance) != 13:
                not_rf.append(folder + '/' + file)
                continue
    
            num_json = num_json + 1
            strainID = resistance[0][0].split('/')[-1]
            test.append(strainID)
                            
            lineage = 'unfilled'
            
            # Find .vcf file if there is one
            vcf_files = list(filter(lambda x: (strainID in x and '.vcf' in x), files))
            
            # Find lineage.txt file if there is one
            lineage_files = list(filter(lambda x: (strainID in x and 'lineage.txt' in x), files))
            
            # If there is a single vcf file, add this to the list to be later converted to lineages
            if len(vcf_files) == 1:
                vcf_dict[strainID] = '//n//groups//gentb_www//predictData//' + folder + '//' + vcf_files[0]
                lineage = 'See array'
                vcf_check = True
                
            # If there are too many vcf files, add the file to the error list and skip!
            elif len(vcf_files) > 1:
                if folder not in IDs_with_errors['More than one .vcf file, did not check for lineage files'].keys():
                    IDs_with_errors['More than one .vcf file, did not check for lineage files'][folder] = [strainID]
                elif strainID not in IDs_with_errors['More than one .vcf file, did not check for lineage files'].get(folder):
                    IDs_with_errors['More than one .vcf file, did not check for lineage files'][folder].append(strainID)
                continue
            
            # If there are no .vcf files, check if there is existing lineage output.
            else:
                if len(lineage_files) == 1:
                    lineage = '//n//groups//gentb_www//predictData//' + folder + '//' + lineage_files[0]
                    
                # If there are no lineage files or too many lineage files, add the file to the error list and skip
                elif len(lineage_files) == 0:
                    if folder not in IDs_with_errors['No lineage files, no vcf'].keys():
                        IDs_with_errors['No lineage files, no vcf'][folder] = [strainID]
                    elif strainID not in IDs_with_errors['No lineage files, no vcf'].get(folder):
                        IDs_with_errors['No lineage files, no vcf'][folder].append(strainID)
                    continue
                else:
                    if folder not in IDs_with_errors['No vcf, More than one lineage file'].keys():
                        IDs_with_errors['No vcf, More than one lineage file'][folder] = [strainID]
                    elif strainID not in IDs_with_errors['No vcf, More than one lineage file'].get(folder):
                        IDs_with_errors['No vcf, More than one lineage file'][folder].append(strainID)
                    continue
                    
            # Case where the strainID is already present:
            if strainID in target_dict.keys():
                other_lineage_file = target_dict[strainID]['Lineage']
                # If a vcf file has already been found and processed, you're good
                if other_lineage_file == 'See array':
                    continue
                # Our existing file thus probably had a lineage file associated - if our new file has a vcf, replace it and overwrite the old dictionary entry:
                elif vcf_check:
                    vcf_dict[strainID] = '//n//groups//gentb_www//predictData//' + folder + '//' + vcf_files[0]
                # This means there are at least two lineage files associated with this one strainID -- it is possible that a .vcf file will be found later though - for now, skip this file (don't overwrite)
                else:
                    if folder not in IDs_with_errors['This strainID is present in more than one folder'].keys():
                        IDs_with_errors['This strainID is present in more than one folder'][folder] = [strainID]
                    elif strainID not in IDs_with_errors['This strainID is present in more than one folder'].get(folder):
                        IDs_with_errors['This strainID is present in more than one folder'][folder].append(strainID)
                    continue
            else:
                unique_json = unique_json + 1
            
            # Finally, add the entry in if it hasn't been caught by any other filter
            target_dict[strainID] = {'Folder':folder, 'Resistance':file, 'Lineage':lineage}
            
            
            
IDs_with_errors['.json file is empty'] = empty_json
IDs_with_errors['.json is not in random forest format'] = not_rf


# In[ ]:


with open('//home//kin672//Kira_GenTB//Isolating Relevant Files//7-18 Output//relevant_file_names_7-18.json', 'w') as f1:
    json.dump(target_dict, f1,  indent=4)
    
with open('//home//kin672//Kira_GenTB//Isolating Relevant Files//7-18 Output//strainIDs_with_errors_7-18.json', 'w') as f2:
    json.dump(IDs_with_errors, f2,  indent=4)  

vcf_to_lineages = get_lineages(vcf_dict)
vcf_to_lineages.to_csv('//home//kin672//Kira_GenTB//Isolating Relevant Files//7-18 Output//vcf_to_lineages_7-18.csv', index = False)

    
print('Number of unique strainIDs that had nonempty .jsons:' + str(unique_json))
print('Total number of nonempty .jsons submitted:' + str(num_json))




