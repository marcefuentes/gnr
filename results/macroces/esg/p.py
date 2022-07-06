#!/usr/bin/python

filename = 'p'
ftype = 'bubbles' # barsall, barsone, scatter
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['ChooseGrain6', 'a2Seen31', 'wmedian']
top_row = [{'treatment':'p', 'control': 'none'}, {'treatment':'p', 'control': 'none'}, {'treatment':'p', 'control': 'none'}]
bottom_row = [{'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}, {'treatment': 'none', 'control': 'optimal'}]
rows = [top_row]

if ftype == 'barsone': glovalue = {'x': -1.0, 'y': 1.0}
