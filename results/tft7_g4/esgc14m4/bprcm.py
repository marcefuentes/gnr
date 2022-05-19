#!/usr/bin/python

filename = 'bprcm'
ftype = 'barsone' # barsall, barsone, scatter
treatment = 'pr'
control = 'none'
drift = False
movie = False
x = 'ES'
y = 'Given'
z0 = 'ChooseGrain'
z1 = 'MimicGrain'

if ftype == 'barsone':
    x_value = -2.0
    y_value = 1.0

