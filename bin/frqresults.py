#! /usr/bin/env python

from glob import glob
import os
import re
import time

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

# Options

traits = ['a2Seen',
          'ChooseGrain',
          'MimicGrain',
          'w']
titles = ['Effort to get $\it{B}$',
          'Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner',
          'Fitness']
folders = ['given0', 'none', 'p', 'r', 'pr', 'p8r']

movie = False
plotsize = 4

# Add data to figure

def figdata(t, lines):
    for f, folder in enumerate(folders):
        df = dfs[f]
        m = df.Time == t
        df = df.loc[m]
        for r, trait in enumerate(traits):
            for a, alpha in enumerate(alphas):
                for l, loges in enumerate(logess):
                    m = (df.alpha == alpha) & (df.logES == loges)
                    d = df.loc[m]
                    freq_a = [col for col in d.columns if re.match(fr'^{trait}\d+$', col)]
                    y = d.loc[:, freq_a]
                    y = y.values[0]
                    y = y.flatten()
                    lines[f, r, a, l].set_ydata(y)
    if movie:
        fig.texts[2].set_text(f't\n{t}')
    return lines.flatten()

# Get data

def read_file(file):
    df = pd.read_csv(file)
    if not movie:
        df = df.tail(1)
    return df

dfs = np.empty(len(folders), dtype=object) 
for i, folder in enumerate(folders):
    filelist = glob(os.path.join(folder, '*.frq'))
    d = list(map(read_file, filelist))
    dfs[i] = pd.concat(d, ignore_index=True)

df = dfs[1]
ts = df.Time.unique()
t = ts[-1]
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
x = np.arange(64)

# Figure properties

width = plotsize*len(traits)
height = plotsize*len(folders)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
step = int(nr/2)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=len(folders),
                             ncols=len(traits))
axs = np.empty((len(folders),
                len(traits),
                len(alphas),
                len(logess)),
               dtype=object)
lines = np.empty(axs.shape, dtype=object)

for f, folder in enumerate(folders):
    for r, trait in enumerate(traits):
        grid = outergrid[f, r].subgridspec(nrows=nr,
                                           ncols=nc,
                                           wspace=0,
                                           hspace=0)
        axs[f, r] = grid.subplots()

left_x = axs[0, 0, 0, 0].get_position().x0
right_x = axs[-1, -1, -1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0, 0, 0].get_position().y1
bottom_y = axs[-1, -1, -1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.5,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=left_x*0.4,
              y=center_y,
              fontsize=biglabels)

if movie:
    fig.text(right_x,
             bottom_y*0.5,
             f't\n0',
             fontsize=biglabels,
             color='grey',
             ha='right')

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.1)

for f, folder in enumerate(folders):
    for r, trait in enumerate(traits):
        letter = ord('a') + f*len(traits) + r
        axs[f, r, 0, 0].set_title(chr(letter),
                                  fontsize=plotsize*5,
                                  pad = 10,
                                  weight='bold',
                                  loc='left')
        if f == 0:
            axs[f, r, 0, 10].set_title(titles[r],
                                       pad=plotsize*9,
                                       fontsize=plotsize*5)
        for a in range(0, nr, step):
            axs[f, r, a, 0].set(yticks=[0.1], yticklabels=[])
            if r == 0:
                axs[f, r, a, 0].set_ylabel(f'{alphas[a]:.1f}',
                                           rotation='horizontal',
                                           horizontalalignment='right',
                                           verticalalignment='center',
                                           y=0.3,
                                           fontsize=ticklabels)
        for l in range(0, nc, step):
            axs[f, r, -1, l].set(xticks=[32], xticklabels=[])
            if folder == folders[-1]:
                axs[f, r, -1, l].set_xlabel(f'{logess[l]:.0f}',
                                            x=0.3,
                                            fontsize=ticklabels)

# Assign lines to axs

dummy_y = np.zeros_like(x)

for f, folder in enumerate(folders):
    for r, trait in enumerate(traits):
        for a, alpha in enumerate(alphas):
            for l, loges in enumerate(logess):
                lines[f, r, a, l], = axs[f, r, a, l].plot(x, dummy_y)

# Add data and save figure

if movie:
    ani = FuncAnimation(fig, figdata, frames=ts, fargs=(lines,), blit=True)
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    figdata(ts[-1], lines,)
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
