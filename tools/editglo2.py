#! /usr/bin/env python

import os

pattern = '.glo'

# loop through all folders in the current directory
# then loop through all subfolders in each folder
# take the last 3 characters of the subfolder name and store them in float variable `given`
# divide `given` by 100 and store the result in `given`
# then loop through all files matching `pattern` in each subfolder
# open each file and replace the string `Given,0.0` with `Given,given`, where `given` is the value stored in the variable `given` without trailing zeros
# if the folder name contains the letter `p`, replace the string `PartnerChoice,0` with `PartnerChoice,1`
# if the folder name contains the letter `r`, replace the string `Reciprocity,0` with `Reciprocity,1`
# if the folder name contains the number `8`, replace the string `GroupSize,2` with `GroupSize,3`

folders = [f for f in os.listdir(os.getcwd()) if os.path.isdir(f)]
for folder in folders:
    subfolders = [f for f in os.listdir(folder) if os.path.isdir(os.path.join(folder, f))]
    for subfolder in subfolders:
        given = float(subfolder[-3:]) / 100
        for f in os.listdir(os.path.join(folder, subfolder)):
            if f.endswith(pattern):
                with open(os.path.join(folder, subfolder, f), 'r') as file:
                    filedata = file.read()
                    filedata = filedata.replace('Given,0.0', f'Given,{given}')
                    if 'p' in folder:
                        filedata = filedata.replace('PartnerChoice,0', 'PartnerChoice,1')
                    if 'r' in folder:
                        filedata = filedata.replace('Reciprocity,0', 'Reciprocity,1')
                    if '8' in folder:
                        filedata = filedata.replace('GroupSize,2', 'GroupSize,3')
                with open(os.path.join(folder, subfolder, f), 'w') as file:
                    file.write(filedata)
