#!/usr/bin/python

filename = 'rcmTL'
ftype = 'barsone' # barsall, barsone, scatter
treatment = 'r'
control = 'none'
drift = False
movie = False
x = 'MimicCost'
y = 'ChooseCost'
z0 = 'ChooseGrain'
z1 = 'MimicGrain'

if ftype == 'barsone':
    x_value = -14.0
    y_value = -7.0
