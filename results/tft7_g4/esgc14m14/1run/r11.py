#!/usr/bin/python

filename = 'r11'
ftype = 'barsone' # barsall, barsone, bubbles, scatterall, scatterone
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['MimicGrain', 'help', 'w']
folderss = [{'treatment':'r', 'control': 'none'}, {'treatment': 'none', 'control': 'r'}]

if ftype == 'barsone': glovalue = {'x': 1.0, 'y': 1.0}
