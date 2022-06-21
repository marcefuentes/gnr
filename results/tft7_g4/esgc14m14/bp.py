#!/usr/bin/python

filename = 'bp0'
ftype = 'barsone' # barsall, barsone, scatter
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['a2Seenmedian', 'helpmedian', 'wmedian']
folderss = [{'treatment':'optimal', 'control': 'none'}, {'treatment': 'none', 'control': 'optimal'}]

if ftype == 'barsone': glovalue = {'x': 0.0, 'y': 1.0}
