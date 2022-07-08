#!/usr/bin/python

filename = 'slope'
ftype = 'bubbles' # barsall, barsone, scatter
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['wmedian', 'wmedian']
top_row = [{'treatment':'optimal', 'control': 'none'}, {'treatment':'optimal', 'control': 'none'}]
bottom_row = [{'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}]
rows = [top_row]

if ftype == 'barsone': glovalue = {'x': -1.0, 'y': 1.0}
