#!/usr/bin/python

filename = 'ir'
ftype = 'scatterall' # barsall, barsone, bubbles, scatterall, scatterone
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = [{'x': 'MimicGrain', 'y': 'w'}, {'x': 'a2Default', 'y': 'w'}, {'x': 'a2Default', 'y': 'MimicGrain'}]
folders = ['r', 'none']
zcolor = 'a2Seen'

if ftype == 'scatterone': glovalue = {'x': 4.0, 'y': 0.9}
