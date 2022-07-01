#!/usr/bin/python

filename = 'pyr'
ftype = 'bubbles' # barsall, barsone, scatter
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['ChooseGrainmedian', 'wmedian', 'wmedian']
top_traits = ['ChooseGrainmedian', 'wmedian', 'wmedian']
bottom_traits = ['MimicGrainmedian', 'wmedian', 'wmedian']
top_row = [{'treatment':'p', 'control': 'none'}, {'treatment':'p', 'control': 'none'}, {'treatment':'p', 'control': 'optimal'}]
bottom_row = [{'treatment': 'r', 'control': 'none'}, {'treatment': 'r', 'control': 'none'}, {'treatment': 'r', 'control': 'optimal'}]
rows = [top_row, bottom_row]

if ftype == 'barsone': glovalue = {'x': -1.0, 'y': 1.0}
