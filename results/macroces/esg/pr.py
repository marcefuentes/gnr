#!/usr/bin/python

filename = 'pr'
ftype = 'bubbles' # barsall, barsone, scatter
movie = False

glos = {'x': 'alpha', 'y': 'Given'}
traits = ['ChooseGrainmedian', 'MimicGrainmedian', 'wmedian']
top_row = [{'treatment':'pr', 'control': 'none'}, {'treatment':'pr', 'control': 'none'}, {'treatment':'pr', 'control': 'optimal'}]
bottom_row = [{'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}]
rows = [top_row]

if ftype == 'barsone': glovalue = {'x': -1.0, 'y': 1.0}
