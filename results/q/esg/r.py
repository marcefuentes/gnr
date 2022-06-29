#!/usr/bin/python

filename = 'r'
ftype = 'bubbles' # barsall, barsone, scatter
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['MimicGrainmedian', 'wmedian', 'wmedian']
top_row = [{'treatment':'r', 'control': 'none'}, {'treatment':'r', 'control': 'none'}, {'treatment':'r', 'control': 'optimal'}]
bottom_row = [{'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}]
rows = [top_row]

if ftype == 'barsone': glovalue = {'x': -1.0, 'y': 1.0}
