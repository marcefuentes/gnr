#!/usr/bin/python

filename = 'r'
ftype = 'bubbles' # barsall, barsone, scatter
movie = False

glos = {'x': 'alpha', 'y': 'Given'}
traits = ['MimicGrain6', 'a2Seen51', 'wmedian']
top_row = [{'treatment':'r', 'control': 'none'}, {'treatment':'r', 'control': 'none'}, {'treatment':'r', 'control': 'none'}]
bottom_row = [{'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}]
rows = [top_row]

if ftype == 'barsone': glovalue = {'x': -1.0, 'y': 1.0}
