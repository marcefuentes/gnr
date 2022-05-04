#!/usr/bin/python

ftype = 'scatter'
drift = False
movie = True
x = 0
y = 1
z0 = 0
z1 = 1
filename = 'scatter'

# Options

xy_list = ('ES', 'Given')
label_list = ('Substitutability of resources', '\nProportion of resource $\it{A}$\ngiven to partner')
log_list = (True, False)
z_list = ('ChooseGrainmedian', 'MimicGrainmedian')
title_list = ('Sensitivity for comparing\npotential partners', 'Sensitivity for comparing\npartner to self')

# Internal

x_axis = xy_list[x]
y_axis = xy_list[y]
x_label = label_list[x]
y_label = label_list[y]
logx = log_list[x]
logy = log_list[y]
c_name_roots = (z_list[z0], z_list[z1])
titles = (title_list[z0], title_list[z1])

# Configuration

x_min = None
x_max = None
y_min = None
y_max = None
bubble_size = 600.0
width = 12.0
height = 6.2
