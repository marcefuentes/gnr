#!/usr/bin/python

filename = 'b0'
ftype = 'barsone' # barsall, barsone, scatter
movie = True
x = 'ES'
y = 'Given'

zg = {}
zg['name'] = 'a2Seen'
zg['treatment'] = 'p'
zg['control'] = 'none'
zw = {}
zw['name'] = 'w'
zw['treatment'] = 'p'
zw['control'] = 'optimal'

z = [zg, zw]

dirs = ['none', 'optimal', 'p']

if ftype == 'barsone':
    x_value = 0.0
    y_value = 1.0
