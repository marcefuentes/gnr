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
titles = ['Prisoner\'s dilemma',
          'Snowdrift',
          'Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner']
traits = ['ChooseGrainmean',
          'MimicGrainmean']
vmaxs = [my.aBmax,
              my.aBmax]
folders = ['given100', 'given095', 'given050']
subfolders = ['p', 'r']

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

df = dfss[0][0]
ts = df.Time.unique()
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1.0 - 1.0/pow(2.0, logess)
x = np.linspace(0.0, my.aBmax, num=ext)
y = np.flip(x)
X, Y = np.meshgrid(x, y)

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

axsss = []
axds = []
for g, row in enumerate(rows):
    axss = []
    for p in range(2):
        grid = outergrid[g, p].subgridspec(nrows=nr,
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
                axs[i, j].set(xticklabels=[], yticklabels=[])
                for axis in ['top','bottom','left','right']:
                    axs[i, j].spines[axis].set_linewidth(0.1)

        for i in range(0, nr, step):
            axs[i, 0].set(yticks=[ext/2]) 
            if p == 0:
                axs[i, 0].set_yticklabels([f'{alphas[i]:3.1f}'],
                                          rotation='horizontal',
                                          horizontalalignment='right',
                                          verticalalignment='center',
                                          y=0.3,
                                          fontsize=ticklabels)
        for j in range(0, nc, step):
            axs[-1, j].set(xticks=[ext/2]) 
            if g == 2:
                axs[-1, j].set_xticklabels([f'{logess[j]:2.0f}'],
                                           x=0.0,
                                           fontsize=ticklabels)
        if g == 0:
            axs[0, 10].set_title(titles[p],
                                 fontsize=plotsize*5,
                                 pad=plotsize*9,
                                 loc='center')
        axss.append(axs)
    axsss.append(axss)

    axd = []
    for c in range(2):
        grid = outergrid[g, c+2].subgridspec(nrows=1,
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
            ax.set_title(titles[c+2],
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

            RR = np.full([ext, ext], rho)
            T = my.fitness(Y, X, given, AA, RR)
            R = my.fitness(Y, Y, given, AA, RR)
            P = my.fitness(X, X, given, AA, RR)
            S = my.fitness(X, Y, given, AA, RR)

            Z = np.full([*T.shape, 4], my.colormap['white'])
            Z = my.prisonercolors(T, R, P, S, Z)
            Z[X >= Y] = [0.9, 0.9, 0.9, 1.0]
            axsss[g][0][i][j].imshow(Z, extent=extent)

            if g == 0:
                Z = np.full([*T.shape, 4], my.colormap['white'])
                Z = my.greycolors(T, R, P, S, Z)
                Z[X >= Y] = [0.9, 0.9, 0.9, 1.0]
                axsss[g][1][i][j].imshow(Z, extent=extent)
            else:
                Z = np.full([*T.shape, 4], my.colormap['white'])
                Z = my.snowdriftcolors(T, R, P, S, Z)
                Z[X >= Y] = [0.9, 0.9, 0.9, 1.0]
                axsss[g][1][i][j].imshow(Z, extent=extent)

    for j, trait in enumerate(traits):
        df = dfss[g][j]
        m = df.Time == ts[-1]
        df = df.loc[m]
        Z = pd.pivot_table(df,
                           values=trait,
                           index='alpha',
                           columns='logES')
        Z = Z.sort_index(axis=0, ascending=False)
        axds[g][j].imshow(Z, vmin=0, vmax=vmaxs[j])

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
