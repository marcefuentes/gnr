#! /usr/bin/env python

from glob import glob
import os
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

# Options

traits = ['ChooseGrainmean',
          'MimicGrainmean']
titles = ['Games',
          'Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner']
vmaxs = [my.a2max,
         my.a2max]
folders = ['given100', 'given95', 'given50']
subfolders = ['none', 'p', 'r']

ext = 256
plotsize = 6

# Data

dfs = np.empty((len(folders), len(subfolders)), dtype=object)
for i, folder in enumerate(folders):
    for j, subfolder in enumerate(subfolders):
        filelist = glob(os.path.join(folder, subfolder, '*.csv'))
        dfs[i, j] = pd.concat(map(pd.read_csv, filelist),
                             ignore_index=True)
        for trait in traits:
            dfs[i, j][trait] = 1.0 - dfs[i, j][trait]

filelist = glob(os.path.join('given00', 'none', '*.csv'))
dfsocial = pd.concat(map(pd.read_csv, filelist),
                     ignore_index=True)

df = dfs[0, 0]
ts = df.Time.unique()
m = dfsocial.Time == ts[-1]
dfso = dfsocial.loc[m]
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1.0 - 1.0/pow(2.0, logess)

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
extent = 0, ext, 0, ext
width = plotsize*len(titles)
height = plotsize*len(folders)
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
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

outergrid = fig.add_gridspec(nrows=len(folders),
                             ncols=len(titles))

axs = np.empty((len(folders),
                len(alphas),
                len(rhos)),
                dtype=object)
axd = np.empty((len(folders),
               len(traits)), 
               dtype=object)

letter = ord('a')
letterposition = 1.035
for g, folder in enumerate(folders):
    grid = outergrid[g, 0].subgridspec(nrows=nr,
                                       ncols=nc,
                                       wspace=0,
                                       hspace=0)
    axs[g] = grid.subplots()
    axs[g, 0, 0].set_title(chr(letter),
                           fontsize=plotsize*5,
                           pad = 10,
                           weight='bold',
                           loc='left')
    letter += 1
    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            axs[g, i, j].set(xticks=[], yticks=[])
            for axis in ['top','bottom','left','right']:
                axs[g, i, j].spines[axis].set_linewidth(0.1)
    if g == 0:
        axs[g, 0, 10].set_title(titles[0],
                                pad=plotsize*9,
                                fontsize=plotsize*5)
    for i in range(0, nr, step):
        axs[g, i, 0].set(yticks=[ext/2]) 
        axs[g, i, 0].set_yticklabels([f'{alphas[i]:3.1f}'],
                                      rotation='horizontal',
                                      horizontalalignment='right',
                                      verticalalignment='center',
                                      y=0.3,
                                      fontsize=ticklabels)
    for j in range(0, nc, step):
        axs[g, -1, j].set(xticks=[ext/2], xticklabels=[]) 
        if g == 2:
            axs[g, -1, j].set_xticklabels([f'{logess[j]:2.0f}'],
                                       x=0.0,
                                       fontsize=ticklabels)

    for c in range(2):
        grid = outergrid[g, c+1].subgridspec(nrows=1,
                                             ncols=1)
        axd[g, c] = grid.subplots()
        axd[g, c].text(0,
                       letterposition,
                       chr(letter),
                       transform=axd[g, c].transAxes,
                       fontsize=plotsize*5,
                       weight='bold')
        letter += 1
        axd[g, c].set(xticks=xticks, yticks=yticks)
        axd[g, c].set(xticklabels=[], yticklabels=[])
        if g == 0:
            axd[g, c].set_title(titles[c+1],
                         pad=plotsize*9,
                         fontsize=plotsize*5)
        if g == 2:
            axd[g, c].set_xticklabels(xticklabels,
                                      x=0.0,
                                      fontsize=ticklabels)

for g, folder in enumerate(folders):

    df = dfs[g][0]
    m = df.Time == ts[-1]
    df = df.loc[m]
    given = df.Given.iloc[0]

    for i, alpha in enumerate(alphas):
        AA = np.full([ext, ext], alpha)
        for j, (rho, loges) in enumerate(zip(rhos, logess)):

            #m = (df.alpha == alpha) & (df.logES == loges)
            #xmin = df.loc[m].a2Seenmean
            #m = (dfso.alpha == alpha) & (dfso.logES == loges)
            #xmax = dfso.loc[m].a2Seenmean
            xmin = 0.0
            xmax = my.a2max
            ymin = xmin
            ymax = xmax
            x = np.linspace(xmin, xmax, num=ext)
            y = np.linspace(ymax, ymin, num=ext)
            X, Y = np.meshgrid(x, y)
            RR = np.full([ext, ext], rho)
            T = my.fitness(Y, X, given, AA, RR)
            R = my.fitness(Y, Y, given, AA, RR)
            P = my.fitness(X, X, given, AA, RR)
            S = my.fitness(X, Y, given, AA, RR)

            Z = my.gamecolors(T, R, P, S)
            Z[X >= Y] = [0.9, 0.9, 0.9, 1.0]
            axs[g, i, j].imshow(Z, extent=extent)

    for j, trait in enumerate(traits):
        df = dfs[g, j + 1]
        m = df.Time == ts[-1]
        df = df.loc[m]
        Z = pd.pivot_table(df,
                           values=trait,
                           index='alpha',
                           columns=['logES'])
        Z = Z.sort_index(axis=0, ascending=False)
        axd[g, j].imshow(Z, vmin=0, vmax=vmaxs[j])

# Save figure

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
