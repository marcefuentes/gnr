#!/usr/bin/python

filename = 'barsall'
ftype = 'barsall' # barsall, barsone, scatter
treatment = 'r'
control = 'none'
movie = False
x = 'ES'
y = 'Given'
z0 = {}
z0['name'] = 'help'
z0['treatment'] = 'r'
z0['control'] = 'none'
z1 = {}
z1['name'] = 'w'
z1['treatment'] = 'r'
z1['control'] = 'none'
z = [z0, z1]

dirs = ['r', 'none']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
