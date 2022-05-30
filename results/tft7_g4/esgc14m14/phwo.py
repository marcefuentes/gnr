#!/usr/bin/python

filename = 'phwo'
ftype = 'scatter' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'

zh = {}
zh['name'] = 'helpmedian'
zh['treatment'] = 'p'
zh['control'] = 'optimal'
zw = {}
zw['name'] = 'wmedian'
zw['treatment'] = 'p'
zw['control'] = 'optimal'

z = [zh, zw]

dirs = ['optimal', 'p']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
