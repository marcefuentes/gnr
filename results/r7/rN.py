#!/usr/bin/python

filename = 'rN'
ftype = 'scatter' # barsall, barsone, scatter
movie = False
x = 'N'
y = 'MimicCost'

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

z02 = {}
z02['name'] = z0['name']
z02['treatment'] = 'none'
z02['control'] = 'optimal'

z12 = {}
z12['name'] = z1['name']
z12['treatment'] = 'none'
z12['control'] = 'optimal'

z22 = {}
z22['name'] = z2['name']
z22['treatment'] = 'none'
z22['control'] = 'optimal'

z = [[z0, z1, z2]]

dirs = ['none', 'optimal', 'r']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0

