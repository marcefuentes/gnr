#!/usr/bin/python

import glob
import imageio
import os
import pandas as pd
import matplotlib.pyplot as plt

# Definitions

x_axis = 'MimicGrainInit'
y_axis = 'ChooseGrainInit'
x_label = 'Sensitivity for changing the level of help'
y_label = 'Sensitivity for switching to a new partner'
c_name_roots = ['a2Seen', 'w']
titles = ['Level of help', 'w']
dirs = ['', 'drift/']

# Properties

fs = 16 # Label font size
sc = 600.0 # Largest bubble
dc = 0.0 # Darkest color
lc = 1.0 # Lightest color

# Data

dfs = []

for d in dirs:
    df = pd.concat(map(pd.read_csv, glob.glob(os.path.join(d, '*.csv'))), ignore_index=True)
    dfs.append(df)

BINS = int(sum(map(lambda x: c_name_roots[0] in x, [*dfs[0]]))/2) - 1

def calculate_colors(z):
    z = lc - (lc - dc)*z
    z = ['%.3f' % elem for elem in z]
    return z

outfile=os.path.splitext(os.path.basename(__file__))[0] + '.gif'
filenames = []

for t in dfs[0].Time.unique():

    print(f'Processing time {t}', end='\r')

    # Figure

    fig, axs = plt.subplots(nrows=len(dfs), ncols=len(c_name_roots), sharex=True, sharey=True, figsize=(16, 10))
    fig.supxlabel(x_label, fontsize=fs)
    fig.supylabel(y_label, fontsize=fs)
    fig.text(0.89, 0.95, f'Time = {t}', fontsize=fs, ha='right')

    # Plots

    for ax, title in zip(axs[0], titles):
        ax.set_title(title, fontsize=fs)

    for row, df in zip(axs, dfs):
        for ax, c_name in zip(row, c_name_roots):
            ax.set_aspect('equal', 'box')
            for b in range(BINS):
                cn = c_name + str(BINS - b - 1)
                ax.scatter(x=df.loc[df.Time == t, x_axis], y=df.loc[df.Time == t, y_axis], c=calculate_colors(df.loc[df.Time == t, cn]), s=(1.00-b/BINS)*(1.00-b/BINS)*sc)

    # Save and close

    filename = f'{outfile}{t}.png'
    filenames.append(filename)

    plt.savefig(filename, dpi=100)
    plt.close()

# Build gif

with imageio.get_writer(outfile, mode='I') as writer:
    for filename in filenames:
        print(f'Adding file {filename} to movie', end='\r')
        image = imageio.imread(filename)
        writer.append_data(image)

for filename in set(filenames):
    os.remove(filename)
