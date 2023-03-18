#! /usr/bin/env python

from glob import glob
import os
import re
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

givens = [1.0, 0.95, 0.5]
titles = ['Effort\nto get B',
          'Fitness']
traits = ['a2Seen',
          'w']
folders = ['given100', 'given95', 'given50']
subfolders = ['none', 'p', 'r']

rows = folders
plotsize = 6

dfss = []
for folder in folders:
    dfs = []
    for subfolder in subfolders:
        filelist = glob(os.path.join(folder, subfolder, '*.frq'))
        df = pd.concat(map(pd.read_csv, filelist),
                       ignore_index=True)
        dfs.append(df)
    dfss.append(dfs)

filelist = glob(os.path.join('given00', 'none', '*.frq'))

df = dfss[0][0]
ts = df.Time.unique()
alphas = np.sort(pd.unique(df.alpha))[::-1]
rowindex = 'alpha'
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
x = np.arange(64)

step = int(nr/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
letterposition = 1.035
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
width = plotsize*len(titles)
height = plotsize*len(rows)
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.515,
              y=0.04,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.04,
              y=0.495,
              fontsize=biglabels)

outergrid = fig.add_gridspec(nrows=len(rows),
                             ncols=len(titles))
axsss = []
for g, row in enumerate(rows):
    axss = []
    for c, trait in enumerate(traits):
        grid = outergrid[g, c].subgridspec(nrows=nr,
                                           ncols=nc,
                                           wspace=0,
                                           hspace=0)
        axs = grid.subplots()
        axs[0, 0].set_title(chr(letter),
                            fontsize=plotsize*5,
                            pad = 10,
                            weight='bold',
                            loc='left')
        letter += 1
        for i, alpha in enumerate(alphas):
            for j, loges in enumerate(logess):
                axs[i, j].set(xticks=[], yticks=[])
                axs[i, j].set(xticklabels=[], yticklabels=[])
                for axis in ['top','bottom','left','right']:
                    axs[i, j].spines[axis].set_linewidth(0.1)
        if c == 0:
            for i in range(0, nr, step):
                axs[i, 0].set(yticks=[0.5]) 
                axs[i, 0].set_yticklabels([f'{alphas[i]:3.1f}'],
                                     rotation='horizontal',
                                     horizontalalignment='right',
                                     verticalalignment='center',
                                     y=0.3,
                                     fontsize=ticklabels)

        if g == 0:
            axs[0, 10].set_title(titles[c],
                                fontsize=plotsize*5,
                                pad=plotsize*9,
                                loc='center')
        if row == rows[-1]:
            for j in range(0, nc, step):
                axs[-1, j].set(xticks=[32], xticklabels=[]) 
                if g == 2:
                    axs[-1, j].set_xticklabels([f'{logess[j]:2.0f}'],
                                               x=0.0,
                                               fontsize=ticklabels)
        axss.append(axs)
    axsss.append(axss)

for g, given in enumerate(givens):

    df = dfss[g][0]
    m = df.Time == ts[-1]
    df = df.loc[m]
    given = df.Given.iloc[0]

    for r, trait in enumerate(traits):
        for i, alpha in enumerate(alphas):
            for j, loges in enumerate(logess):
                m = (df.alpha == alpha) & (df.logES == loges)
                d = df.loc[m]
                freq_a = [col for col in d.columns if re.match(fr'^{trait}\d+$', col)]
                y = d.loc[:, freq_a]
                y = y.values[0]
                y = y.flatten()
                axsss[g][r][i][j].plot(x, y)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
