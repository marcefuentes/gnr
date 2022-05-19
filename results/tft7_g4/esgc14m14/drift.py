#!/usr/bin/python

filename = 'drift'
ftype = 'scatter' # barsall, barsone, scatter
treatment = 'pr'
control = 'none'
drift = True
movie = False
x = 'ES'
y = 'Given'
z0 = 'helpmedian'
z1 = 'wmedian'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
