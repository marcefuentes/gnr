#!/usr/bin/python

filename = 'brcm'
ftype = 'barsone' # barsall, barsone, scatter
treatment = 'r'
control = 'none'
drift = False
movie = False
z0 = 'ChooseGrain'
z1 = 'MimicGrain'

if ftype == 'barsone':
    x_value = -15.0
    y_value = -5.0
else:
    x = 'N'
    y = 'ChooseCost'

