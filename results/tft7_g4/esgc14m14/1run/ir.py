#!/usr/bin/python

filename = 'ir_MimicGrain_insensitive_-51Seen'
sensitive = False
ftype = 'scatterone' # barsall, barsone, bubbles, scatterall, scatterone
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = [{'x': 'MimicGrain', 'y': 'w'}, {'x': 'a2Seen', 'y': 'w'}, {'x': 'MimicGrain', 'y': 'a2Seen'}]
folders = ['r', 'none']
zalpha = 'MimicGrain'

if ftype == 'scatterone': glovalue = {'x': -5.0, 'y': 1.0}
