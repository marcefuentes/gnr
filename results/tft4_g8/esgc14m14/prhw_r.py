#!/usr/bin/python

filename = 'prhw_r'
ftype = 'scatter' # barsall, barsone, scatter
treatment = 'pr'
control = 'p'
drift = False
movie = False
x = 'ES'
y = 'Given'
z0 = 'helpmedian'
z1 = 'wmedian'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0