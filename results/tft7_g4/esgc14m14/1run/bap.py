#!/usr/bin/python

filename = 'bap'
ftype = 'barsall' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'

z0 = {}
z0['name'] = 'ChooseGrain'
z0['treatment'] = 'p'
z0['control'] = 'none'

z1 = {}
z1['name'] = 'a2Seen'
z1['treatment'] = 'p'
z1['control'] = 'none'

z2 = {}
z2['name'] = 'w'
z2['treatment'] = 'p'
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

dirs = ['none', 'optimal', 'p']

if (ftype == 'barsone') or (ftype == 'barsall'):
    ymax = 0.1

if ftype == 'barsone':
    x_value = -5.0
    y_value = 1.0

