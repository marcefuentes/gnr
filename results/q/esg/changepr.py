#!/usr/bin/python

filename = 'changepr'
ftype = 'scatter' # barsall, barsone, scatter
movie = False
x = 'alpha'
y = 'Given'

z0 = {}
z0['name'] = 'chose_partner'
z0['treatment'] = 'pr'
z0['control'] = 'none'

z1 = {}
z1['name'] = 'changed_a2'
z1['treatment'] = 'pr'
z1['control'] = 'none'

z2 = {}
z2['name'] = 'wmedian'
z2['treatment'] = 'pr'
z2['control'] = 'optimal'

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

z = [[z0, z1]]

dirs = ['none', 'pr']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0

