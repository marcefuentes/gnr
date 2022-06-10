#!/usr/bin/python

filename = 'br-5'
ftype = 'barsone' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'

z0 = {}
z0['name'] = 'MimicGrain'
z0['treatment'] = 'r'
z0['control'] = 'none'

z1 = {}
z1['name'] = 'MimicGrain'
z1['treatment'] = 'none'
z1['control'] = 'optimal'

z2 = {}
z2['name'] = 'w'
z2['treatment'] = 'r'
z2['control'] = 'none'

z02 = {}
z02['name'] = 'MimicGrainmedian'
z02['treatment'] = 'r'
z02['control'] = 'none'

z12 = {}
z12['name'] = 'wmedian'
z12['treatment'] = 'r'
z12['control'] = 'none'

z22 = {}
z22['name'] = z2['name']
z22['treatment'] = 'r'
z22['control'] = 'optimal'

z = [[z0, z1, z2]]

dirs = ['none', 'optimal', 'r']

if ftype == 'barsone':
    x_value = -5.0
    y_value = 1.0

