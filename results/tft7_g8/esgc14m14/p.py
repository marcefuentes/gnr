#!/usr/bin/python

filename = 'p'
ftype = 'scatter' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'

zg = {}
zg['name'] = 'ChooseGrainmedian'
zg['treatment'] = 'p'
zg['control'] = 'none'
zw = {}
zw['name'] = 'wmedian'
zw['treatment'] = 'p'
zw['control'] = 'optimal'

z = [zg, zw]

dirs = ['none', 'optimal', 'p']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
