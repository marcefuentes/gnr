#!/usr/bin/python

filename = 'r'
ftype = 'scatter' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'

zg = {}
zg['name'] = 'MimicGrainmedian'
zg['treatment'] = 'r'
zg['control'] = 'none'
zw = {}
zw['name'] = 'wmedian'
zw['treatment'] = 'r'
zw['control'] = 'optimal'

z = [zg, zw]

dirs = ['none', 'optimal', 'r']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
