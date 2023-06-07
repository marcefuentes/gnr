#! /usr/bin/env python

import os

givens = ['given000', 'given050', 'given095', 'given100']
lmhs = ['low', 'mid', 'high']
nprs = ['none', 'p', 'p8', 'p8r', 'pr', 'r']

# Create a list of npr names in the current npr that start with discrete and have a number after
discretes = []
others = []

for f in os.listdir('.'):
    if f.startswith('discrete') and f[8].isdigit():
        discretes.append(f)
    else:
        others.append(f)

for discrete in discretes:
    for i, lmh in enumerate(lmhs):
        os.mkdir(discrete + str(i))

for discrete in discretes:
    for given in givens:
        if os.path.exists(discrete + '/' + given):
            for i, lmh in enumerate(lmhs):
                if os.path.exists(discrete + '/' + given + '/' + lmh):
                    if not os.path.exists(discrete + str(i)):
                        os.mkdir(discrete + str(i))
                    if not os.path.exists(discrete + str(i) + '/' + given):
                        os.mkdir(discrete + str(i) + '/' + given)
                    for npr in nprs:
                        if os.path.exists(discrete + '/' + given + '/' + lmh + '/' + npr):
                            if not os.path.exists(discrete + str(i) + '/' + npr):
                                os.mkdir(discrete + str(i) + '/' + npr)
                            if not os.path.exists(discrete + str(i) + '/' + npr + '/' + given):
                                os.mkdir(discrete + str(i) + '/' + npr + '/' + given)
                            if os.listdir(discrete + '/' + given + '/' + lmh + '/' + npr):
                                os.system('mv ' + discrete + '/' + given + '/' + lmh + '/' + npr + '/*.* ' + discrete + str(i) + '/' + npr + '/' + given + '/')
                            os.system('rmdir ' + discrete + '/' + given + '/' + lmh + '/' + npr)
                    os.system('rmdir ' + discrete + '/' + given + '/' + lmh)
            os.system('rmdir ' + discrete + '/' + given)
    os.system('rmdir ' + discrete)
                      
for other in others:        
    for npr in nprs:
        os.mkdir(other + '/' + npr)
        for given in givens:
            os.mkdir(other + '/' + npr + '/' + given)

    # Move files to the correct nprs
    for given in givens:
        if os.path.exists(other + '/' + given):
            for npr in nprs:
                if os.path.exists(other + '/' + given + '/' + npr):
                    if os.listdir(other + '/' + given + '/' + npr):
                        os.system('mv ' + other + '/' + given + '/' + npr + '/*.* ' + other + '/' + npr + '/' + given + '/')
                    os.system('rmdir ' + other + '/' + given + '/' + npr)
            os.system('rmdir ' + other + '/' + given)
    
