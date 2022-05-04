#!/usr/bin/python

ftype = 'barsone' # barsall, barsone, scatter
drift = False
movie = False
x = 0
y = 1
z0 = 9
z1 = 10
x_value = -1.0
y_value = 1.0
filename = 'scatter'

# Options

xy_list = []
label_list = []
log_list = []

# xy 0
xy_list.append('ES')
label_list.append('Substitutability of resources')
log_list.append(True)
# xy 1
xy_list.append('Given')
label_list.append('\nProportion of resource $\it{A}$\ngiven to partner')
log_list.append(False)
# xy 2
xy_list.append('ChooseCost')
label_list.append('Cost of comparing potential partners')
log_list.append(True)
# xy 3
xy_list.append('MimicCost')
label_list.append('Cost of comparing partner to self')
log_list.append(True)
# xy 4
xy_list.append('MimicCost')
label_list.append('Cost of comparing potential partners to self')
log_list.append(True)
# xy 5
xy_list.append('DeathRate')
label_list.append('Death rate')
log_list.append(True)
# xy 6
xy_list.append('GroupSize')
label_list.append('Potential partners')
log_list.append(True)
# xy 7
xy_list.append('MimicGrainInit')
label_list.append('Sensitivity for comparing partner to self')
log_list.append(False)
# xy 8
xy_list.append('MimicGrainInit')
label_list.append('Sensitivity for comparing potential partners to self')
log_list.append(False)
# xy 9
xy_list.append('ChooseGrainInit')
label_list.append('Sensitivity for comparing potential partners')
log_list.append(False)

z_list = []
title_list = []

# z 0
z_list.append('ChooseGrainmedian')
title_list.append('Sensitivity for comparing\npotential partners')
# z 1
z_list.append('MimicGrainmedian')
title_list.append('Sensitivity for comparing\npartner to self')
# z 2
z_list.append('MimicGrainmedian')
title_list.append('Sensitivity for comparing\npotential partners to self')
# z 3
z_list.append('Helpmedian')
title_list.append('Help')
# z 4
z_list.append('wmedian')
title_list.append('Fitness')
# z 5
z_list.append('chose_partner')
title_list.append('Frequency of\nswitching to a new partner')
# z 6
z_list.append('changed_a2')
title_list.append('Frequency of\nchanging help')
# z 7
z_list.append('HelpBD')
title_list.append('Fluctuation of help')
# z 8
z_list.append('wBD')
title_list.append('Fluctuation of fitness')
# z 9
z_list.append('ChooseGrain')
title_list.append('Sensitivity for comparing\npotential partners')
# z 10
z_list.append('MimicGrain')
title_list.append('Sensitivity for comparing\npartner to self')
# z 11
z_list.append('MimicGrain')
title_list.append('Sensitivity for comparing\npotential partners to self')

# Internal

x_axis = xy_list[x]
x_label = label_list[x]
logx = log_list[x]
y_axis = xy_list[y]
y_label = label_list[y]
logy = log_list[y]
c_name_roots = (z_list[z0], z_list[z1])
titles = (title_list[z0], title_list[z1])

# Configuration

x_min = None
x_max = None
y_min = None
y_max = None
if xy_list[x] == 'DeathRate':
    x_min = 0.005
    x_max = 0.2
if xy_list[y] == 'GroupSize':
    y_min = 48.0
    y_max = 2.6
if 'bars' in ftype:
    y_max = 0.2
bubble_size = 600.0
if (('chose' in z_list[z0]) or ('BD' in z_list[z0])):
    bubble_size = 2000.0
width = 12.0
height = 6.2
