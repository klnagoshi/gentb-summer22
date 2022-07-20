# Testing new method for pulling strain IDs -- is it better than the last method?
import pandas as pd
import numpy as np
import json
import os

result = {}
empty_json = []
predictData = os.listdir('//n//groups//gentb_www//predictData')

for folder in predictData:
    files = os.listdir('//n//groups//gentb_www//predictData//' + folder)
    
    for file in files:
        if 'matrix.json' in file:
            # PULL STRAIN ID FROM FIRST ENTRY IN EACH OUTPUT JSON
            strainID1 = file.replace('.matrix.json','') # OLD METHOD
            json_path = '//n//groups//gentb_www//predictData//' + folder + '//' + file
            if open(json_path).read() == '':
                empty_json.append(folder + '/' + file)
            else:
                with open(json_path) as f:
                    resistance = json.load(f)
                strainID2 = resistance[0][0][0].split('/')[-1]

                if strainID1 != strainID2:
                    if folder in result.keys():
                        result[folder].append([strainID1, strainID2])
                    else:
                        result[folder]=([strainID1, strainID2])
print(empty_json)                
print(json.dumps(result))