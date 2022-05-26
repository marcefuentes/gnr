#!/usr/bin/python

filename = 'pw'
ftype = 'scatter' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'
z0 = {}
z0['name'] = 'ChooseGrainmedian'
z0['treatment'] = 'p'
z0['control'] = 'none'
z1 = {}
z1['name'] = 'wmedian'
z1['treatment'] = 'p'
z1['control'] = 'none'
z2 = {}
z2['name'] = 'wmedian'
z2['treatment'] = 'p'
z2['control'] = 'optimal'
z = [z0, z1, z2]

dirs = ['p', 'none', 'optimal']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
