#!/usr/bin/python

filename = 'optimal'
ftype = 'scatter' # barsall, barsone, scatter
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['a2Seenmedian', 'helpmedian', 'wmedian']
folderss = [{'treatment':'optimal', 'control': 'none'}, {'treatment': 'none', 'control': 'optimal'}]

if ftype == 'barsone': glovalue = {'x': -1.0, 'y': 1.0}
