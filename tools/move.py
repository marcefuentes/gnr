#! /usr/bin/env python

import os

folders = ['none', 'p', 'p8', 'p8r', 'pr', 'r']
subfolders = ['given000', 'given050', 'given095', 'given100']

# Create folders and the 4 subfolders in each
for folder in folders:
    os.mkdir(folder)
    for subfolder in subfolders:
        os.mkdir(folder + '/' + subfolder)

# Move files to the correct folders
for subfolder in subfolders:
    for folder in folders:
        # Move all files *.* to folder/subfolder/
        os.system('mv ' + subfolder + '/' + folder + '/*.* ' + folder + '/' + subfolder + '/')
        # Remove folder/subfolder/
        os.system('rmdir ' + subfolder + '/' + folder)
    # Remove subfolder/
    os.system('rmdir ' + subfolder)
    
