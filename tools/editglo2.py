#! /usr/bin/env python

import os

pattern = '.glo'
given_folder = os.path.basename(os.getcwd())
decimal_part = given_folder[-3:]
given = float(decimal_part) / 100
folders = [f for f in os.listdir(os.getcwd()) if os.path.isdir(f)]

for folder in folders:
    for f in os.listdir(folder):
        if f.endswith(pattern):
            with open(os.path.join(folder, f), 'r') as file:
                filedata = file.read()
                filedata = filedata.replace('Given,0.0', f'Given,{given}')
                if 'p' in folder:
                    filedata = filedata.replace('PartnerChoice,0', 'PartnerChoice,1')
                if 'r' in folder:
                    filedata = filedata.replace('Reciprocity,0', 'Reciprocity,1')
                if '8' in folder:
                    filedata = filedata.replace('GroupSize,2', 'GroupSize,3')
            with open(os.path.join(folder, f), 'w') as file:
                file.write(filedata)
