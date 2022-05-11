#!/usr/bin/python

filename = 'rhw'
ftype = 'scatter' # barsall, barsone, scatter
treatment = 'r'
control = 'none'
drift = False
movie = False
x = 'MimicCost'
y = 'ChooseCost'
z0 = 'helpmedian'
z1 = 'wmedian'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
