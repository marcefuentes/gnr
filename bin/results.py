#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

traits = ['a2Seenmean',
          'ChooseGrainmean',
          'MimicGrainmean',
          'wmean']
titles = ['Effort to get $\it{B}$',
          'Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner',
          'Fitness']
traitvmaxs = [my.a2max, my.a2max, my.a2max, 2.0]
folders = ['given0', 'none', 'p', 'r', 'pr', 'p8r']

movie = False
rows = folders
plotsize = 4

def adddata(t, ims):
    for i, df in enumerate(dfs):
        for j, trait in enumerate(traits):
            m = df.Time == t
            df = df.loc[m]
            Z = pd.pivot_table(df,
                               values=trait,
                               index=[rowindex],
                               columns=['logES'])
            Z = Z.sort_index(axis=0, ascending=False)
            ims[i, j].set_array(Z) 
    if movie:
        fig.texts[2].set_text(f't\n{t}')
    return ims.flatten()

dfs = []
for folder in folders:
    filelist = glob(os.path.join(folder, '*.csv'))
    df = pd.concat(map(pd.read_csv, filelist),
                   ignore_index=True)
    for trait in traits:
        if 'Grain' in trait:
            df[trait] = 1.0 - df[trait]
    dfs.append(df)

df = dfs[1]
ts = df.Time.unique()
if movie:
    frames = []
else:
    ts = [ts[-1]]
rowindex = 'alpha'
nr = df['alpha'].nunique()
nc = df['logES'].nunique()

xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
letterposition = 1.035
xticks = [0, nc/2-0.5, nc-1]
yticks = [0, nr/2-0.5, nr-1]
xmin = df['logES'].min()
xmax = df['logES'].max()
ymin = df['alpha'].min()
ymax = df['alpha'].max()
xticklabels = [f'{xmin:2.0f}',
               f'{(xmin + xmax)/2.0:2.0f}',
               f'{xmax:2.0f}']
yticklabels = [f'{ymax:3.1f}',
               f'{(ymin + ymax)/2.0:3.1f}',
               f'{ymin:3.1f}']
width = plotsize*len(titles)
height = plotsize*len(rows)
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig, axs = plt.subplots(nrows=len(rows),
                        ncols=len(traits),
                        figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.515,
              y=0.06,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.03,
              y=0.493,
              fontsize=biglabels)
if movie:
    fig.text(0.90,
             0.93,
             f't\n0',
             fontsize=biglabels,
             color='grey',
             ha='right')
for ax in fig.get_axes():
    ax.set(xticks=xticks, yticks=yticks)
    ax.set(xticklabels=[], yticklabels=[])
    ax.text(0,
            letterposition,
            chr(letter),
            transform=ax.transAxes,
            fontsize=plotsize*5,
            weight='bold')
    letter += 1
for i, row in enumerate(rows):
    axs[i, 0].set_yticklabels(yticklabels, fontsize=ticklabels)
for j, title in enumerate(titles):
    axs[0, j].set_title(title, pad=plotsize*10, fontsize=plotsize*5)
    axs[-1, j].set_xticklabels(xticklabels, fontsize=ticklabels)

ims = np.empty(axs.shape, dtype=object) 
dummy_Z = np.empty((nr, nc), dtype=np.float32)

for i, df in enumerate(dfs):
    for j, trait in enumerate(traits):
        ims[i, j] = axs[i, j].imshow(dummy_Z,
                                     vmin=0,
                                     vmax=traitvmaxs[j])
if movie:
    ani = FuncAnimation(fig, adddata, frames=ts, fargs=(ims,), blit=True)
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    adddata(ts[-1], ims,)
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
