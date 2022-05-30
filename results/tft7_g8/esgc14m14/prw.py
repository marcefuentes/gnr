#!/usr/bin/python

filename = 'prw'
ftype = 'scatter' # barsall, barsone, scatter
movie = True
x = 'ES'
y = 'Given'

zc = {}
zc['name'] = 'ChooseGrainmedian'
zc['treatment'] = 'pr'
zc['control'] = 'none'
zm = {}
zm['name'] = 'MimicGrainmedian'
zm['treatment'] = 'pr'
zm['control'] = 'none'
zwn = {}
zwn['name'] = 'wmedian'
zwn['treatment'] = 'pr'
zwn['control'] = 'none'
zwo = {}
zwo['name'] = 'wmedian'
zwo['treatment'] = 'pr'
zwo['control'] = 'optimal'

z = [zc, zm, zwn, zwo]

dirs = ['none', 'optimal', 'pr']

if ftype == 'barsone':
    x_value = -1.0
    y_value = 1.0
