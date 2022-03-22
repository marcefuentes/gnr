#!/usr/bin/python

import glob
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

# Definitions

p_names = ['ChooseCost', 'MimicCost']
c_name_roots = ['ChooseGrain', 'MimicGrain', 'a2Seen']
x_label = 'Time'
y_label = 'Frequency'
titles = ['Sensitivity for\nswitching to a new partner', 'Sensitivity for\nchanging the level of help', 'Level of help']
dirs = ['', 'drift/']

# Properties

fs = 16 # Label font size


# Data

if len(sys.argv) < 3:
    print(f'You must enter two arguments: {p_names[0]} {p_names[1]}')
    exit(1)

p_values = [float(str("{:.6f}".format(pow(2, -int(sys.argv[1]))))), float(str("{:.6f}".format(pow(2, -int(sys.argv[2])))))]

outfile=os.path.splitext(os.path.basename(__file__))[0] + '_' + sys.argv[1] + '_' + sys.argv[2] + '.png'

dfs = []

for d in dirs:
    df = pd.concat(map(pd.read_csv, glob.glob(os.path.join(d, '*.csv'))), ignore_index=True)
    df = df[(df[p_names[0]] == p_values[0]) & (df[p_names[1]] == p_values[1])]
    dfs.append(df)

BINS = int(sum(map(lambda x: c_name_roots[0] in x, [*dfs[0]]))/2) - 1

color_list_g = []
color_list_h = []

for b in range(BINS):
    c = b/BINS
    color_list_g.append('%.3f' % c)
    color_list_h.append('%.3f' % (1.0 - c))

color_lists = [color_list_g, color_list_g, color_list_h]

# Figure

fig, axs = plt.subplots(nrows=len(dfs), ncols=len(c_name_roots), sharex=True, sharey=True, figsize=(16, 10))
fig.supxlabel(x_label, fontsize=fs)
fig.supylabel(y_label, fontsize=fs)

# Plots

for ax, title in zip(axs[0], titles):
    ax.set_title(title, fontsize=fs)

for row, df in zip(axs, dfs):
    for ax, c_name, color_list in zip(row, c_name_roots, color_lists):
        y_data = []
        for b in range(BINS):
            y_data.append(df[c_name + str(b)])
        ax.stackplot(df.Time, y_data, colors=color_list)

# Save and close

plt.savefig(outfile)
plt.close()

