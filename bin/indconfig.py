#!/usr/bin/python

filename = 'ir_MimicGrain_sensitive'
ftype = 'scatterall' # barsall, barsone, bubbles, scatterall, scatterone
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = [{'x': 'MimicGrain', 'y': 'w'}, {'x': 'a2Default', 'y': 'w'}, {'x': 'a2Default', 'y': 'MimicGrain'}]
folders = ['r', 'none']
zalpha = 'MimicGrain'
sensitive = True

if ftype == 'scatterone': glovalue = {'x': 0.0, 'y': 1.0}
