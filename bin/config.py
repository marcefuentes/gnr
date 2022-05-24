#!/usr/bin/python

filename = 'grain'
ftype = 'scatter' # barsall, barsone, scatter
drift = False
movie = False
x = 'ES'
y = 'Given'
z0 = 'ChooseGrainmedian'
treatment0 = 'pr'
control0 = 'none'
z1 = 'MimicGrainmedian'
treatment1 = 'pr'
control1= 'none'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
