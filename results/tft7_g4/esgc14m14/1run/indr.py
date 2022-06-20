#!/usr/bin/python

filename = 'bard'
ftype = 'barsall' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'

z00 = {}
z00['x'] = 'a2Seen'
z00['y'] = 'w'
z00['treatment'] = 'r'

z01 = {}
z01['x'] = 'MimicGrain'
z01['y'] = 'w'
z01['treatment'] = 'r'

z02 = {}
z02['x'] = 'a2Seen'
z02['y'] = 'MimicGrain'
z02['treatment'] = 'r'

z10 = {}
z10['x'] = z00['x']
z10['y'] = z00['y']
z10['treatment'] = 'none'

z11 = {}
z11['x'] = z01['x']
z11['y'] = z01['y']
z11['treatment'] = 'none'

z12 = {}
z12['x'] = z02['x']
z12['y'] = z02['y']
z12['treatment'] = 'none'

z = [[z00, z01, z02], [z10, z11, z12]]

dirs = ['none', 'r']

if (ftype == 'barsone') or (ftype == 'barsall'):
    ymax = 0.1

if ftype == 'barsone':
    x_value = -5.0
    y_value = 1.0

