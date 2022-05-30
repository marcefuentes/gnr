#!/usr/bin/python

filename = 'pr'
ftype = 'scatter' # barsall, barsone, scatter
movie = False
x = 'ES'
y = 'Given'

zcg = {}
zcg['name'] = 'ChooseGrainmedian'
zcg['treatment'] = 'pr'
zcg['control'] = 'none'
zrg = {}
zrg['name'] = 'MimicGrainmedian'
zrg['treatment'] = 'pr'
zrg['control'] = 'none'
zw = {}
zw['name'] = 'wmedian'
zw['treatment'] = 'pr'
zw['control'] = 'optimal'

z = [zcg, zrg, zw]

dirs = ['none', 'optimal', 'pr']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
