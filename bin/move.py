#! /usr/bin/env python

import os

folders = ['given000', 'given050', 'given095', 'given100']
new_starting_names = ['000', '441', '882', '1323']
subfolders = ['none', 'p', 'p8', 'p8r', 'pr', 'r']
extensions = ['csv', 'glo', 'gl2', 'frq']
# Each folder has all subfolders
# For example, 'given000' has 'none', 'p', 'p8', 'p8r', 'pr', 'r' subfolders
# Each subfolder has files with 4 extensions: csv, glo, gl2, frq
# For each extension there are 441 files, named 101.* - 541.*

# Create folder 'data' if it does not exist
if not os.path.exists('data'):
    os.makedirs('data')
# Create subfolders in 'data'
for subfolder in subfolders:
    if not os.path.exists('data/' + subfolder):
        os.makedirs('data/' + subfolder)

for folder, new_starting_name in zip(folders, new_starting_names):
    for subfolders in subfolders:
        # Rename files in subfolder starting with new_starting_names
        # For example, 'given000/none/101.csv' is renamed to 'data/none/000.csv'
        # and 'given050/none/101.csv' is renamed to 'data/none/441.csv'
        for extension in extensions:
            for i in range(101, 542):
                old_name = folder + '/' + subfolder + '/' + str(i) + '.' + extension
                new_name = 'data/' + subfolder + '/' + str(i - 100 + int(new_starting_name)).zfill(3) + '.' + extension

