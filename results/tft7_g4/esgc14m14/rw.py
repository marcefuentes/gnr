#!/usr/bin/python

filename = 'rw'
ftype = 'scatter' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'
z0 = {}
z0['name'] = 'MimicGrainmedian'
z0['treatment'] = 'r'
z0['control'] = 'none'
z1 = {}
z1['name'] = 'wmedian'
z1['treatment'] = 'r'
z1['control'] = 'none'
z2 = {}
z2['name'] = 'wmedian'
z2['treatment'] = 'r'
z2['control'] = 'optimal'
z = [z0, z1, z2]

dirs = ['r', 'none', 'optimal']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
