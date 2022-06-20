#!/usr/bin/python

filename = 'optimal'
ftype = 'scatter' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'

z0 = {}
z0['name'] = 'a2Seenmedian'
z0['treatment'] = 'none'
z0['control'] = 'optimal'

z1 = {}
z1['name'] = 'helpmedian'
z1['treatment'] = 'none'
z1['control'] = 'optimal'

z2 = {}
z2['name'] = 'wmedian'
z2['treatment'] = 'none'
z2['control'] = 'optimal'

z02 = {}
z02['name'] = z0['name']
z02['treatment'] = 'optimal'
z02['control'] = 'none'

z12 = {}
z12['name'] = z1['name']
z12['treatment'] = 'optimal'
z12['control'] = 'none'

z22 = {}
z22['name'] = z2['name']
z22['treatment'] = 'optimal'
z22['control'] = 'none'

z = [[z0, z1, z2], [z02, z12, z22]]

dirs = ['none', 'optimal']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0

