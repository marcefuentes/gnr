#!/usr/bin/python

filename = 'rcm'
ftype = 'scatter' # barsall, barsone, scatter
treatment = 'r'
control = 'none'
drift = False
movie = False
x = 'ES'
y = 'Given'
z0 = 'ChooseGrainmedian'
z1 = 'MimicGrainmedian'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
