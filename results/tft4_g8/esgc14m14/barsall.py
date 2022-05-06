#!/usr/bin/python

filename = 'grainbarsall'
ftype = 'barsall' # barsall, barsone, scatter
treatment = 'pr'
control = 'none'
drift = False
movie = False
x = 'ES'
y = 'Given'
z0 = 'ChooseGrain'
z1 = 'MimicGrain'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
