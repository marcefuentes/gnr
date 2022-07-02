#!/usr/bin/python

filename = 'bar'
ftype = 'barsall' # barsall, barsone, bubbles, scatterall, scatterone
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['MimicGrain', 'a2Seen', 'w']
top_row = [{'treatment':'r', 'control': 'none'}, {'treatment': 'r', 'control': 'none'}, {'treatment': 'r', 'control': 'none'}]
rows = [top_row]

if ftype == 'barsone': glovalue = {'x': -1.0, 'y': 1.0}
