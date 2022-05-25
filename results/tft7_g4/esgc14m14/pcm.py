#!/usr/bin/python

filename = 'pcm'
ftype = 'scatter' # barsall, barsone, scatter
drift = False
movie = False
x = 'ES'
y = 'Given'
z0 = {}
z0['name'] = 'ChooseGrainmedian'
z0['treatment'] = 'p'
z0['control'] = 'none'
z1['name'] = 'MimicGrainmedian'
z1['treatment'] = 'p'
z1['control'] = 'none'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
