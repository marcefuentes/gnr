#!/usr/bin/python

filename = 'lossnone'
ftype = 'scatter' # barsall, barsone, scatter
treatment = 'none'
control = 'optimal'
drift = False
movie = False
x = 'ES'
y = 'Given'
z0 = 'helpmedian'
z1 = 'wmedian'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0