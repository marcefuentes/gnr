#!/usr/bin/python

filename = 'r2'
ftype = 'bubbles' # barsall, barsone, scatter
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['MimicGrain6', 'a2Seen31', 'wmedian']
top_row = [{'treatment':'r', 'control': 'none'}, {'treatment':'r', 'control': 'none'}, {'treatment':'r', 'control': 'none'}]
bottom_row = [{'treatment': 'r', 'control': 'optimal'}, {'treatment': 'r', 'control': 'optimal'}, {'treatment': 'r', 'control': 'optimal'}]
rows = [top_row, bottom_row]

if ftype == 'barsone': glovalue = {'x': -1.0, 'y': 1.0}
