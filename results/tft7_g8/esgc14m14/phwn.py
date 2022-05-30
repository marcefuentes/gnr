#!/usr/bin/python

filename = 'phwn'
ftype = 'scatter' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'

zh = {}
zh['name'] = 'helpmedian'
zh['treatment'] = 'p'
zh['control'] = 'none'
zw = {}
zw['name'] = 'wmedian'
zw['treatment'] = 'p'
zw['control'] = 'none'

z = [zh, zw]

dirs = ['none', 'p']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
