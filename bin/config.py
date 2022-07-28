#!/usr/bin/python

filename = 'pyr'
ftype = 'bubbles' # barsall, barsone, scatter
movie = False

glos = {'x': 'ES', 'y': 'Given'}
top_traits = ['ChooseGrainmedian', 'wmedian', 'wmedian']
bottom_traits = ['MimicGrainmedian', 'wmedian', 'wmedian']
traits = [top_traits, bottom_traits]
top_folders = [{'treatment':'p', 'control': 'none'}, {'treatment':'p', 'control': 'none'}, {'treatment':'p', 'control': 'optimal'}]
bottom_folders = [{'treatment': 'r', 'control': 'none'}, {'treatment': 'r', 'control': 'none'}, {'treatment': 'r', 'control': 'optimal'}]
folders = [top_folders, bottom_folders]

if ftype == 'barsone':
    glovalue = {'x': -5.0, 'y': 0.9}
    ylim = 1.0

if ftype == 'barsall':
    ylim = 1.0
