#! /usr/bin/env python

import os
import fileinput
import sys

pattern = '*.glo'
abs_folder_path = os.path.abspath(os.getcwd())
parent_folder = os.path.basename(os.path.dirname(abs_folder_path))
folder_name = os.path.basename(abs_folder_path)
decimal_part = parent_folder[-3:]
given = float(decimal_part) / 100

for line in fileinput.input(files=[os.path.join(abs_folder_path, f) for f in os.listdir(abs_folder_path) if f.endswith('.glo')], inplace=True, backup='',):
    sys.stdout.write(line.replace('Given,0.0', f'Given,{given}'))

if 'p' in folder_name:
    for line in fileinput.input(files=[os.path.join(abs_folder_path, f) for f in os.listdir(abs_folder_path) if f.endswith('.glo')], inplace=True, backup='',):
        sys.stdout.write(line.replace('PartnerChoice,0', 'PartnerChoice,1'))

if 'r' in folder_name:
    for line in fileinput.input(files=[os.path.join(abs_folder_path, f) for f in os.listdir(abs_folder_path) if f.endswith('.glo')], inplace=True, backup='',):
        sys.stdout.write(line.replace('Reciprocity,0', 'Reciprocity,1'))

if '8' in folder_name:
    for line in fileinput.input(files=[os.path.join(abs_folder_path, f) for f in os.listdir(abs_folder_path) if f.endswith('.glo')], inplace=True, backup='',):
        sys.stdout.write(line.replace('GroupSize,2', 'Groupsize,3'))

