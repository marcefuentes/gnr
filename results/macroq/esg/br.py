#!/usr/bin/python

filename = 'br50'
ftype = 'barsone' # barsall, barsone, scatter
movie = False

glos = {'x': 'alpha', 'y': 'Given'}
traits = ['MimicGrain', 'a2Seen', 'w']
top_row = [{'treatment':'r', 'control': 'none'}, {'treatment':'r', 'control': 'none'}, {'treatment':'r', 'control': 'none'}]
bottom_row = [{'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}]
rows = [top_row, bottom_row]

if ftype == 'barsone': glovalue = {'x': 5.0, 'y': 0.0}
