#!/usr/bin/python

filename = 'indr0-5'
ftype = 'scatterone' # barsall, barsone, bubbles, scatterall, scatterone
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = [{'x': 'MimicGrain', 'y': 'w'}, {'x': 'a2Seen', 'y': 'w'}, {'x': 'a2Seen', 'y': 'MimicGrain'}]
folders = ['r', 'none']

if ftype == 'scatterone': glovalue = {'x': -5.0, 'y': 0.0}
