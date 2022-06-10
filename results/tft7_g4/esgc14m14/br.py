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
z0['max'] = 0.1

z1 = {}
z1['name'] = 'a2Seen'
z1['treatment'] = 'r'
z1['control'] = 'none'
z1['max'] = 0.1

z2 = {}
z2['name'] = 'w'
z2['treatment'] = 'r'
z2['control'] = 'none'
z2['max'] = 0.1

z02 = {}
z02['name'] = z0['name']
z02['treatment'] = 'none'
z02['control'] = 'optimal'
z02['max'] = z0['max']

z12 = {}
z12['name'] = z1['name']
z12['treatment'] = 'none'
z12['control'] = 'optimal'
z12['max'] = z1['max']

z22 = {}
z22['name'] = z2['name']
z22['treatment'] = 'none'
z22['control'] = 'optimal'
z22['max'] = z2['max']

z = [[z0, z1, z2], [z02, z12, z22]]

dirs = ['none', 'optimal', 'r']

if ftype == 'barsone':
    x_value = -5.0
    y_value = 1.0

