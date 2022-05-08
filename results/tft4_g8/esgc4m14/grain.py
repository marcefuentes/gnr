#!/usr/bin/python

filename = 'grain'
ftype = 'scatter' # barsall, barsone, scatter
drift = False
movie = False
x = 'ES'
y = 'Given'
z0 = 'ChooseGrainmedian'
z1 = 'MimicGrainmedian'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
