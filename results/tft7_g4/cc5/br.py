#!/usr/bin/python

filename = 'br14'
ftype = 'barsone' # barsall, barsone, scatter
movie = True
x = 'MimicCost'
y = 'ChooseCost'

z0 = {}
z0['name'] = 'MimicGrain'
z0['treatment'] = 'r'
z0['control'] = 'none'

z1 = {}
z1['name'] = 'MimicGrain'
z1['treatment'] = 'none'
z1['control'] = 'r'

z2 = {}
z2['name'] = 'w'
z2['treatment'] = 'r'
z2['control'] = 'none'

z02 = {}
z02['name'] = z0['name']
z02['treatment'] = 'none'
z02['control'] = 'r'

z12 = {}
z12['name'] = z1['name']
z12['treatment'] = 'none'
z12['control'] = 'r'

z22 = {}
z22['name'] = z2['name']
z22['treatment'] = 'optimal'
z22['control'] = 'r'

z = [[z0, z1, z2]]

dirs = ['none', 'r']

if ftype == 'barsone':
    x_value = -14.0
    y_value = -14.0

