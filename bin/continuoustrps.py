#! /usr/bin/env python

from glob import glob
from matplotlib import cm
import matplotlib.pyplot as plt
import mymodule as my
import numpy as np
import os
import pandas as pd
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

trait = 'ChooseGrainmean'
folders = ['given100', 'given95', 'given50']
subfolders = ['none', 'p']

rows = folders
plotsize = 8

dfss = []
for folder in folders:
    dfs = []
    for subfolder in subfolders:
        filelist = glob(os.path.join(folder, subfolder, '*.csv'))
        df = pd.concat(map(pd.read_csv, filelist),
                        ignore_index=True)
        df[trait] = 1.0 - df[trait]
        dfs.append(df)
    dfss.append(dfs)

filelist = glob(os.path.join('given00', 'none', '*.csv'))
dfsocial = pd.concat(map(pd.read_csv, filelist),
                     ignore_index=True)

df = dfss[0][0]
ts = df.Time.unique()
alphas = np.sort(pd.unique(df.alpha))[::-1]
rowindex = 'alpha'
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1.0 - 1.0/pow(2.0, logess)
RR, AA = np.meshgrid(rhos, alphas)
df = dfsocial.loc[df.Time == ts[-1]]
highs = pd.pivot_table(df,
                       values='a2Seenmean',
                       index=[rowindex],
                       columns=['logES'])
highs = highs.sort_index(axis=0, ascending=False)
highs = highs.to_numpy()

xlim=[0, 5]
ylim=[0.0, 2.0]
step = int(nr/2)
xaxis = [1, 2, 3, 4]
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
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
width = plotsize
height = plotsize*len(rows)
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(width*1.25, height))
fig.supxlabel(xlabel,
              x=0.502,
              y=0.04,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.04,
              y=0.5,
              fontsize=biglabels)

outergrid = fig.add_gridspec(nrows=len(rows),
                             ncols=1,
                             left=0.22,
                             right=0.80)

axss = []
for g, row in enumerate(rows):
    grid = outergrid[g].subgridspec(nrows=nr,
                                    ncols=nc,
                                    wspace=0,
                                    hspace=0)
    axs = grid.subplots()
    axs[0, 0].text(0,
                   4.8,
                   chr(letter),
                   fontsize=plotsize*5,
                   weight='bold')
    letter += 1

    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            axs[i, j].set(xticks=[], yticks=[])
            axs[i, j].set(xlim=xlim, ylim=ylim)
            for axis in ['top','bottom','left','right']:
                axs[i, j].spines[axis].set_linewidth(0.1)
    for i in range(0, nr, step):
        axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                             rotation='horizontal',
                             horizontalalignment='right',
                             verticalalignment='center',
                             fontsize=ticklabels)
    if g == 2:
        for j in range(0, nc, step):
            axs[-1, j].set_xlabel(f'{logess[j]:2.0f}',
                                  x=0.3,
                                  fontsize=ticklabels)
    axss.append(axs)

for g, folder in enumerate(folders):

    df = dfss[g][0]
    df = df.loc[df.Time == ts[-1]]
    given = df.Given.iloc[0]
    lows = pd.pivot_table(df,
                          values='a2Seenmean',
                          index=[rowindex],
                          columns=['logES'])
    lows = lows.sort_index(axis=0, ascending=False)
    lows = lows.to_numpy()
    T = my.fitness(highs, lows, given, AA, RR)
    R = my.fitness(highs, highs, given, AA, RR)
    P = my.fitness(lows, lows, given, AA, RR)
    S = my.fitness(lows, highs, given, AA, RR)

    df = dfss[g][1]
    df = df.loc[df.Time == ts[-1]]
    Z = pd.pivot_table(df,
                       values=trait,
                       index=[rowindex],
                       columns=['logES'])
    Z = Z.sort_index(axis=0, ascending=False)
    Z = Z.to_numpy()

    axs = axss[g]
    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            y = [T[i, j], R[i, j], P[i, j], S[i, j]]
            for line in axs[i, j].get_lines():
                line.remove()
            axs[i, j].plot(xaxis,
                           y,
                           c=cm.viridis(Z[i, j]),
                           linewidth=1,
                           marker='o',
                           markerfacecolor='white',
                           markersize=plotsize/3)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
