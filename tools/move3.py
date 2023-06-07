#! /usr/bin/env python

import os

folders = ['none', 'p', 'p8', 'p8r', 'pr', 'r']
subfolders = ['given000', 'given050', 'given095', 'given100']
discretefolders = ['low', 'mid', 'high']
for f in os.listdir('.'):
    if f.startswith('discrete') and f[8].isdigit():
        for subfolder in subfolders:
            if os.path.exists(f + '/' + subfolder):
                for i, discretefolder in enumerate(discretefolders):
                    if os.path.exists(f + '/' + subfolder + '/' + discretefolder):
                        if not os.path.exists(f + str(i)):
                            os.mkdir(f + str(i))
                        for folder in folders:
                            if os.path.exists(f + '/' + subfolder + '/' + discretefolder + '/' + folder):
                                if not os.path.exists(f + str(i) + '/' + folder):
                                    os.mkdir(f + str(i) + '/' + folder)
                                    if not os.path.exists(f + str(i) + '/' + folder + '/' + subfolder):
                                        os.mkdir(f + str(i) + '/' + folder + '/' + subfolder)
                          

for f in os.listdir('.'):
    if f.startswith('discrete') and f[8].isdigit():
        for subfolder in subfolders:
            if os.path.exists(f + '/' + subfolder):
                for i, discretefolder in enumerate(discretefolders):
                    for folder in folders:
                        if os.path.exists(f + '/' + subfolder + '/' + discretefolder + '/' + folder):
                            if os.listdir(f + '/' + subfolder + '/' + discretefolder + '/' + folder):
                                os.system('mv ' + f + '/' + subfolder + '/' + discretefolder + '/' + folder + '/*.* ' + f + str(i) + '/' + folder + '/' + subfolder + '/')
                            os.system('rmdir ' + f + '/' + subfolder + '/' + discretefolder + '/' + folder)
                    os.system('rmdir ' + f + '/' + subfolder + '/' + discretefolder)
                os.system('rmdir ' + f + '/' + subfolder)
        
    else:
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
    
