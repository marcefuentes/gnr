#!/usr/bin/python

filename = 'ip_ChooseGrain_insensitive'
sensitive = False
ftype = 'scatterall' # barsall, barsone, bubbles, scatterall, scatterone
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = [{'x': 'ChooseGrain', 'y': 'w'}, {'x': 'a2Default', 'y': 'w'}, {'x': 'a2Default', 'y': 'ChooseGrain'}]
top_row = 'p'
bottom_row = 'none'
rows = [top_row, bottom_row]
zalpha = 'ChooseGrain'

if ftype == 'scatterone': glovalue = {'x': -5.0, 'y': 1.0}
