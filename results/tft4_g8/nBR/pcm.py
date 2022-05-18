#!/usr/bin/python

filename = 'pcm'
ftype = 'scatter' # barsall, barsone, scatter
treatment = 'p'
control = 'none'
drift = False
movie = False
x = 'N'
y = 'ChooseCost'
z0 = 'ChooseGrainmedian'
z1 = 'MimicGrainmedian'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
