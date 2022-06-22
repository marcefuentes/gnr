#!/usr/bin/python

filename = 'indr'
ftype = 'scatterall' # barsall, barsone, bubbles, scatterall, scatterone
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = [{'x': 'MimicGrain', 'y': 'w'}, {'x': 'a2Seen', 'y': 'w'}, {'x': 'a2Seen', 'y': 'MimicGrain'}]
folders = ['r', 'none']

if ftype == 'scatterone': glovalue = {'x': 0.0, 'y': 1.0}
