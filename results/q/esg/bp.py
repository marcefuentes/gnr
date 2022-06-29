#!/usr/bin/python

filename = 'bp-11'
ftype = 'barsone' # barsall, barsone, scatter
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['ChooseGrain', 'a2Seen', 'w']
top_row = [{'treatment':'p', 'control': 'none'}, {'treatment':'p', 'control': 'none'}, {'treatment':'p', 'control': 'none'}]
bottom_row = [{'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}]
rows = [top_row, bottom_row]

if ftype == 'barsone': glovalue = {'x': -1.0, 'y': 1.0}
