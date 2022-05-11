#!/usr/bin/python

filename = 'praw'
ftype = 'scatter' # barsall, barsone, scatter
treatment = 'pr'
control = 'none'
drift = False
movie = False
x = 'MimicCost'
y = 'ChooseCost'
z0 = 'a2Seenmedian'
z1 = 'wmedian'

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
