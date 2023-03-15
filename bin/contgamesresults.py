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

givens = [1.0, 0.95, 0.5]
titles = ['Games',
          'Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner']
traits = ['ChooseGrainmean',
          'MimicGrainmean']
traitvmaxs = [my.a2max,
              my.a2max]
folders = ['given100', 'given95', 'given50']
subfolders = ['none', 'p', 'r']

ext = 256
plotsize = 6
rows = folders

dfss = []
for folder in folders:
    dfs = []
    for subfolder in subfolders:
        filelist = glob(os.path.join(folder, subfolder, '*.csv'))
        df = pd.concat(map(pd.read_csv, filelist),
                       ignore_index=True)
        for trait in traits:
            df[trait] = 1.0 - df[trait]
        dfs.append(df)
    dfss.append(dfs)

filelist = glob(os.path.join('given00', 'none', '*.csv'))
dfsocial = pd.concat(map(pd.read_csv, filelist),
                     ignore_index=True)

df = dfss[0][0]
ts = df.Time.unique()
m = dfsocial.Time == ts[-1]
dfso = dfsocial.loc[m]
alphas = np.sort(pd.unique(df.alpha))[::-1]
rowindex = 'alpha'
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1.0 - 1.0/pow(2.0, logess)

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
extent = 0, ext, 0, ext
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

axss = []
axds = []
for g, row in enumerate(rows):
    grid = outergrid[g, 0].subgridspec(nrows=nr,
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
        for j, rho in enumerate(rhos):
            axs[i, j].set(xticks=[], yticks=[])
            for axis in ['top','bottom','left','right']:
                axs[i, j].spines[axis].set_linewidth(0.1)
    for i in range(0, nr, step):
        axs[i, 0].set(yticks=[ext/2]) 
        axs[i, 0].set_yticklabels([f'{alphas[i]:3.1f}'],
                             rotation='horizontal',
                             horizontalalignment='right',
                             verticalalignment='center',
                             y=0.3,
                             fontsize=ticklabels)

    if g == 0:
        axs[0, 10].set_title(titles[0],
                            fontsize=plotsize*5,
                            pad=plotsize*9,
                            loc='center')
    for j in range(0, nc, step):
        axs[-1, j].set(xticks=[ext/2], xticklabels=[]) 
        if g == 2:
            axs[-1, j].set_xticklabels([f'{logess[j]:2.0f}'],
                                       x=0.0,
                                       fontsize=ticklabels)
    axss.append(axs)

    axd = []
    for c in range(2):
        grid = outergrid[g, c+1].subgridspec(nrows=1,
                                             ncols=1)
        ax = grid.subplots()
        ax.text(0,
                letterposition,
                chr(letter),
                transform=ax.transAxes,
                fontsize=plotsize*5,
                weight='bold')
        letter += 1
        ax.set(xticks=xticks, yticks=yticks)
        ax.set(xticklabels=[], yticklabels=[])
        if g == 0:
            ax.set_title(titles[c+1],
                         pad=plotsize*9,
                         fontsize=plotsize*5)
        if g == 2:
            ax.set_xticklabels(xticklabels,
                               x=0.0,
                               fontsize=ticklabels)
        axd.append(ax)
    axds.append(axd)

for g, given in enumerate(givens):

    df = dfss[g][0]
    m = df.Time == ts[-1]
    df = df.loc[m]
    given = df.Given.iloc[0]

    for i, alpha in enumerate(alphas):
        AA = np.full([ext, ext], alpha)
        for j, (rho, loges) in enumerate(zip(rhos, logess)):

            xmin = 0.0
            m = (df.alpha == alpha) & (df.logES == loges)
            xma = df.loc[m].a2Seenmean
            m = (dfso.alpha == alpha) & (dfso.logES == loges)
            ymi = dfso.loc[m].a2Seenmean
            xmax = (ymi - xma)/2.0
            ymin = xmax
            ymax = my.a2max
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
            axss[g][i][j].imshow(Z, extent=extent)

    for j, trait in enumerate(traits):
        df = dfss[g][j + 1]
        m = df.Time == ts[-1]
        df = df.loc[m]
        Z = pd.pivot_table(df,
                           values=trait,
                           index=[rowindex],
                           columns=['logES'])
        Z = Z.sort_index(axis=0, ascending=False)
        axds[g][j].imshow(Z, vmin=0, vmax=traitvmaxs[j])

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
