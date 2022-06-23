#!/usr/bin/python

filename = 'banone'
ftype = 'barsall' # barsall, barsone, scatter
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['a2Seen', 'help', 'w']
folderss = [{'treatment':'none', 'control': 'optimal'}, {'treatment': 'optimal', 'control': 'none'}]

if ftype == 'barsone': glovalue = {'x': 0.0, 'y': 1.0}
