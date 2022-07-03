#!/usr/bin/python

filename = 'ir_MimicGrain_insensitive'
sensitive = False
ftype = 'scatterall' # barsall, barsone, bubbles, scatterall, scatterone
movie = False

glos = {'x': 'alpha', 'y': 'Given'}
traits = [{'x': 'MimicGrain', 'y': 'w'}, {'x': 'a2Default', 'y': 'w'}, {'x': 'a2Default', 'y': 'MimicGrain'}]
top_row = 'r'
bottom_row = 'none'
rows = [top_row, bottom_row]
zalpha = 'MimicGrain'

if ftype == 'scatterone': glovalue = {'x': -5.0, 'y': 1.0}
