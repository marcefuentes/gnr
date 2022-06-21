#!/usr/bin/python

filename = 'bap'
ftype = 'barsall' # barsall, barsone, scatter
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = ['ChooseGrain', 'a2Seen', 'w']
folderss = [{'treatment':'p', 'control': 'none'}, {'treatment': 'none', 'control': 'p'}]

if ftype == 'barsone': glovalue = {'x': 0.0, 'y': 1.0}
