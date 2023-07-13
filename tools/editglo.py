#! /usr/bin/env python

import os

pattern = ".glo"

# loop through all folders ("none", "p", etc.) in the current directory
# then loop through all subfolders ("given000", etc.) in each folder

folders = [f for f in os.listdir(os.getcwd()) if os.path.isdir(f)]
for folder in folders:
    subfolders = [f for f in os.listdir(folder) if os.path.isdir(os.path.join(folder, f))]
    for subfolder in subfolders:
        given = float(subfolder[-3:]) / 100
        for f in os.listdir(os.path.join(folder, subfolder)):
            if f.endswith(pattern):
                with open(os.path.join(folder, subfolder, f), "r") as file:
                    filedata = file.read()
                    filedata = filedata.replace("Given,0.0", f"Given,{given}")
                    if "p" in folder:
                        filedata = filedata.replace("PartnerChoice,0", "PartnerChoice,1")
                    if "r" in folder:
                        filedata = filedata.replace("Reciprocity,0", "Reciprocity,1")
                    if "i" in folder:
                        filedata = filedata.replace("Reciprocity,0", "Reciprocity,1")
                        filedata = filedata.replace("IndirectR,0", "IndirectR,1")
                    if "l" in folder:
                        filedata = filedata.replace("Language,0", "Language,1")
                    if "8" in folder:
                        filedata = filedata.replace("GroupSize,2", "GroupSize,3")
                with open(os.path.join(folder, subfolder, f), "w") as file:
                    file.write(filedata)
