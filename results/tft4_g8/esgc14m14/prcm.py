#!/usr/bin/python

filename = 'prcm'
ftype = 'scatter' # barsall, barsone, scatter
treatment = 'pr'
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
