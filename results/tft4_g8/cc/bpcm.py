#!/usr/bin/python

filename = 'bpcm'
ftype = 'barsall' # barsall, barsone, scatter
treatment = 'p'
control = 'none'
drift = False
movie = False
x = 'MimicCost'
y = 'ChooseCost'
z0 = 'ChooseGrain'
z1 = 'MimicGrain'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
