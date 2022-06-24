#!/usr/bin/python

filename = 'ip_ChooseGrain_sensitive_-51'
ftype = 'scatterone' # barsall, barsone, bubbles, scatterall, scatterone
movie = False

glos = {'x': 'ES', 'y': 'Given'}
traits = [{'x': 'ChooseGrain', 'y': 'w'}, {'x': 'a2Default', 'y': 'w'}, {'x': 'a2Default', 'y': 'ChooseGrain'}]
folders = ['p', 'none']
zalpha = 'ChooseGrain'
sensitive = True

if ftype == 'scatterone': glovalue = {'x': -5.0, 'y': 1.0}
