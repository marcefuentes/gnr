#!/usr/bin/python

filename = 'phw'
ftype = 'scatter' # barsall, barsone, scatter
treatment = 'p'
control = 'none'
drift = False
movie = True
x = 'ES'
y = 'Given'
z0 = 'helpmedian'
z1 = 'wmedian'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
