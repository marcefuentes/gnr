#!/usr/bin/python

filename = 'bp'
ftype = 'barsone' # barsall, barsone, scatter
movie = True
x = 'ES'
y = 'Given'

zc = {}
zc['name'] = 'ChooseGrain'
zc['treatment'] = 'p'
zc['control'] = 'none'
zm = {}
zm['name'] = 'MimicGrain'
zm['treatment'] = 'p'
zm['control'] = 'none'
zwn = {}
zwn['name'] = 'w'
zwn['treatment'] = 'p'
zwn['control'] = 'none'
zwo = {}
zwo['name'] = 'w'
zwo['treatment'] = 'p'
zwo['control'] = 'optimal'

z = [zc, zm, zwn, zwo]

dirs = ['none', 'optimal', 'p']

if ftype == 'barsone':
    x_value = 2.0
    y_value = 1.0
