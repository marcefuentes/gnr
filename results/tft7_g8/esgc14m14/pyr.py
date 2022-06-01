#!/usr/bin/python

filename = 'pyr'
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
z1['control'] = 'optimal'

z2 = {}
z2['name'] = 'wmedian'
z2['treatment'] = 'optimal'
z2['control'] = 'none'

z02 = {}
z02['name'] = 'MimicGrainmedian'
z02['treatment'] = 'r'
z02['control'] = 'none'

z12 = {}
z12['name'] = 'wmedian'
z12['treatment'] = 'r'
z12['control'] = 'optimal'

z22 = {}
z22['name'] = z2['name']
z22['treatment'] = 'none'
z22['control'] = 'optimal'

z = [[z0, z1], [z02, z12]]

dirs = ['none', 'optimal', 'p', 'r']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0

