#!/usr/bin/python

filename = 'bp0'
ftype = 'barsone' # barsall, barsone, scatter
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['MimicGrain', 'a2Seen', 'w']
folderss = [{'treatment':'p', 'control': 'none'}, {'treatment': 'none', 'control': 'p'}]

if ftype == 'barsone': glovalue = {'x': 0.0, 'y': 1.0}
