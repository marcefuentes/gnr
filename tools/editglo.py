#! /usr/bin/env python

import os

input_file_extension = ".glo"

# loop through all mechanisms ("none", "p", etc.) in the current directory
# then loop through all givens ("given000", etc.) in each mechanism

mechanisms = [f for f in os.listdir(os.getcwd()) if os.path.isdir(f)]
for mechanism in mechanisms:
    givens = [f for f in os.listdir(mechanism) if os.path.isdir(os.path.join(mechanism, f))]
    for given in givens:
        given_value = float(given[-3:]) / 100
        for f in os.listdir(os.path.join(mechanism, given)):
            if f.endswith(input_file_extension):
                with open(os.path.join(mechanism, given, f), "r") as file:
                    filedata = file.read()
                    filedata = filedata.replace("Given,0.0", f"Given,{given_value}")
                    if "p" in mechanism:
                        filedata = filedata.replace("PartnerChoice,0", "PartnerChoice,1")
                    if "r" in mechanism:
                        filedata = filedata.replace("Reciprocity,0", "Reciprocity,1")
                    if "i" in mechanism:
                        filedata = filedata.replace("Reciprocity,0", "Reciprocity,1")
                        filedata = filedata.replace("IndirectR,0", "IndirectR,1")
                    if "l" in mechanism:
                        filedata = filedata.replace("Language,0", "Language,1")
                    if "8" in mechanism:
                        filedata = filedata.replace("GroupSize,2", "GroupSize,3")
                with open(os.path.join(mechanism, given, f), "w") as file:
                    file.write(filedata)
