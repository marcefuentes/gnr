#!/usr/bin/python

filename = 'optimal'
ftype = 'scatter' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'

z0 = {}
z0['name'] = 'a2Seenmedian'
z0['treatment'] = 'optimal'
z0['control'] = 'none'

z1 = {}
z1['name'] = 'helpmedian'
z1['treatment'] = 'optimal'
z1['control'] = 'none'

z2 = {}
z2['name'] = 'wmedian'
z2['treatment'] = 'optimal'
z2['control'] = 'none'

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

z = [[z0, z1, z2], [z02, z12, z22]]

dirs = ['none', 'optimal']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0

