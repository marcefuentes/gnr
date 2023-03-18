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
folders = ['given0', 'none']

movie = False
plotsize = 4

# Add data to figure

def figdata(t, lines):
    for i, folder in enumerate(folders):
        df = dfs[i]
        m = df.Time == t
        df = df.loc[m]
        for j, trait in enumerate(traits):
            for k, alpha in enumerate(alphas):
                for l, loges in enumerate(logess):
                    m = (df.alpha == alpha) & (df.logES == loges)
                    d = df.loc[m]
                    freq_a = [col for col in d.columns if re.match(fr'^{trait}\d+$', col)]
                    y = d.loc[:, freq_a]
                    y = y.values[0]
                    y = y.flatten()
                    lines[i, j, k, l].set_ydata(y)
    if movie:
        fig.texts[2].set_text(f't\n{t}')
    return lines.flatten()

# Get data

dfs = np.empty(len(folders), dtype=object) 
for i, folder in enumerate(folders):
    filelist = glob(os.path.join(folder, '*.frq'))
    dfs[i] = pd.concat(map(pd.read_csv, filelist),
                   ignore_index=True)

df = dfs[1]
ts = df.Time.unique()
t = ts[-1]
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
x = np.arange(64)

# Figure properties

step = int(nr/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
xticks = [0, nc/2-0.5, nc-1]
yticks = [0, nr/2-0.5, nr-1]
xmin = logess[0]
xmax = logess[-1]
ymin = alphas[-1]
ymax = alphas[0]
xticklabels = [f'{xmin:2.0f}',
               f'{(xmin + xmax)/2.0:2.0f}',
               f'{xmax:2.0f}']
yticklabels = [f'{ymax:3.1f}',
               f'{(ymin + ymax)/2.0:3.1f}',
               f'{ymin:3.1f}']
width = plotsize*len(traits)
height = plotsize*len(folders)
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.515,
              y=0.06,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.03,
              y=0.493,
              fontsize=biglabels)
outergrid = fig.add_gridspec(nrows=len(folders),
                             ncols=len(traits))
axs = np.empty((len(folders),
                len(traits),
                nr,
                nc),
                dtype=object)
if movie:
    fig.text(0.90,
             0.93,
             f't\n0',
             fontsize=biglabels,
             color='grey',
             ha='right')

letter = ord('a')
letterposition = 4.8
for i, folder in enumerate(folders):
    for j, trait in enumerate(traits):
        grid = outergrid[i, j].subgridspec(nrows=nr,
                                           ncols=nc,
                                           wspace=0,
                                           hspace=0)
        axs[i, j] = grid.subplots()
        axs[i, j, 0, 0].set_title(chr(letter),
                                  fontsize=plotsize*5,
                                  pad = 10,
                                  weight='bold',
                                  loc='left')
        letter += 1
        for k, alpha in enumerate(alphas):
            for l, loges in enumerate(logess):
                axs[i, j, k, l].set(xticks=[], yticks=[])
                for axis in ['top','bottom','left','right']:
                    axs[i, j, k, l].spines[axis].set_linewidth(0.1)
        if i == 0:
            axs[i, j, 0, 10].set_title(titles[j],
                                       pad=plotsize*9,
                                       fontsize=plotsize*5)
        for k in range(0, nr, step):
            axs[i, j, k, 0].set(yticks=[0.1], yticklabels=[])
            if j == 0:
                axs[i, j, k, 0].set_ylabel(f'{alphas[i]:3.1f}',
                                           rotation='horizontal',
                                           horizontalalignment='right',
                                           verticalalignment='center',
                                           y=0.3,
                                           fontsize=ticklabels)
        for l in range(0, nc, step):
            axs[i, j, -1, l].set(xticks=[32], xticklabels=[])
            if folder == folders[-1]:
                axs[i, j, -1, l].set_xlabel(f'{logess[j]:2.0f}',
                                            x=0.3,
                                            fontsize=ticklabels)

lines = np.empty(axs.shape, dtype=object)
dummy_y = np.zeros_like(x)

for i, folder in enumerate(folders):
    for j, trait in enumerate(traits):
        for k, alpha in enumerate(alphas):
            for l, loges in enumerate(logess):
                lines[i, j, k, l], = axs[i, j, k, l].plot(x, dummy_y)

# Save figure

if movie:
    ani = FuncAnimation(fig, figdata, frames=ts, fargs=(lines,), blit=True)
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    figdata(ts[-1], lines,)
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
