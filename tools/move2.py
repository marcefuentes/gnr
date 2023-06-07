#! /usr/bin/env python

import os

folders = ['none', 'p', 'p8', 'p8r', 'pr', 'r']
subfolders = ['given000', 'given050', 'given095', 'given100']

# in all directories in the current directory
for f in os.listdir('.'):
    # Create folders and the 4 subfolders in each
    for folder in folders:
        os.mkdir(f + '/' + folder)
        for subfolder in subfolders:
            os.mkdir(f + '/' + folder + '/' + subfolder)

    # Move files to the correct folders
    for subfolder in subfolders:
        if os.path.exists(f + '/' + subfolder):
            for folder in folders:
                if os.path.exists(f + '/' + subfolder + '/' + folder):
                    if os.listdir(f + '/' + subfolder + '/' + folder):
                        os.system('mv ' + f + '/' + subfolder + '/' + folder + '/*.* ' + f + '/' + folder + '/' + subfolder + '/')
                    os.system('rmdir ' + f + '/' + subfolder + '/' + folder)
            os.system('rmdir ' + f + '/' + subfolder)
    
