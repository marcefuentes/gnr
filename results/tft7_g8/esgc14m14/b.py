#!/usr/bin/python

filename = 'b'
ftype = 'barsone' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'

zg = {}
zg['name'] = 'help'
zg['treatment'] = 'optimal'
zg['control'] = 'none'

zw = {}
zw['name'] = 'w'
zw['treatment'] = 'optimal'
zw['control'] = 'none'

zg2 = {}
zg2['name'] = zg['name']
zg2['treatment'] = 'none'
zg2['control'] = 'optimal'

zw2 = {}
zw2['name'] = zw['name']
zw2['treatment'] = 'none'
zw2['control'] = 'optimal'

z = [[zg, zw], [zg2, zw2]]

dirs = ['none', 'optimal']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 0.0

