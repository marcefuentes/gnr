#!/usr/bin/python

filename = 'r'
ftype = 'bubbles' # barsall, barsone, scatter
movie = True

glos = {'x': 'alpha', 'y': 'Given'}
traits = ['MimicGrain0', 'changed_a2', 'wmedian']
top_row = [{'treatment':'r', 'control': 'none'}, {'treatment':'r', 'control': 'none'}, {'treatment':'r', 'control': 'optimal'}]
bottom_row = [{'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}]
rows = [top_row]

if ftype == 'barsone': glovalue = {'x': -1.0, 'y': 1.0}
