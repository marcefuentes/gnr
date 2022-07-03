#!/usr/bin/python

filename = 'bap'
ftype = 'barsall' # barsall, barsone, bubbles, scatterall, scatterone
movie = False

glos = {'x': 'alpha', 'y': 'Given'}
traits = ['ChooseGrain', 'a2Seen31', 'w']
top_row = [{'treatment':'p', 'control': 'none'}, {'treatment': 'p', 'control': 'none'}, {'treatment': 'p', 'control': 'none'}]
rows = [top_row]

if ftype == 'barsone': glovalue = {'x': -1.0, 'y': 1.0}
