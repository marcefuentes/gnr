#!/usr/bin/python

filename = 'brcm'
ftype = 'barsone' # barsall, barsone, scatter
treatment = 'r'
control = 'none'
drift = False
movie = True
x = 'Given'
y = 'ES'
z0 = 'ChooseGrain'
z1 = 'MimicGrain'

if ftype == 'barsone':
    x_value = 1.0
    y_value = -5.0
