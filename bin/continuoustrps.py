#! /usr/bin/env python

from glob import glob
import matplotlib.pyplot as plt
import mymodule as my
import numpy as np
import os
import pandas as pd
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

titles = ['Games',
          'Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner']
traits = ['ChooseGrainmean',
          'MimicGrainmean']
traitvmaxs = [my.a2max,
              my.a2max]
folders = ['given100', 'given95', 'given50']
subfolders = ['none', 'p', 'r']

rows = folders
plotsize = 4

dfss = []
for folder in folders:
    dfs = []
    for subfolder in subfolders:
        df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, subfolder, '*.csv'))),
                        ignore_index=True)
        df.ChooseGrainmean = 1.0 - df.ChooseGrainmean
        df.MimicGrainmean = 1.0 - df.MimicGrainmean
        dfs.append(df)
    dfss.append(dfs)

dfsocial = pd.concat(map(pd.read_csv, glob(os.path.join('given00', 'none', '*.csv'))),
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
highs = pd.pivot_table(dfsocial.loc[df.Time == ts[-1]],
            values='a2Seenmean',
            index=[rowindex],
            columns=['logES']).sort_index(axis=0,
                                        ascending=False)
highs = highs.to_numpy()
xlim=[0, 5]
ylim=[0.0, 2.0]
step = int(nr/2)
xaxis = [1, 2, 3, 4]
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
ticklabels = plotsize*4
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.513,
              y=0.03,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.04,
              y=0.493,
              fontsize=biglabels)

outergrid = fig.add_gridspec(nrows=len(rows),
                             ncols=len(titles))
                             #left=0.20,
                             #right=0.80)
                             #top=0.8,
                             #bottom=0.2)

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
                        weight='bold',
                        loc='left')
    letter += 1

    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            axs[i, j].set(xticks=[], yticks=[])
            axs[i, j].set(xlim=xlim, ylim=ylim)
    for i in range(0, nr, step):
        axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                             rotation='horizontal',
                             horizontalalignment='right',
                             verticalalignment='center',
                             fontsize=ticklabels)
    if g == 0:
        axs[0, int(nc/2)].set_title(titles[g],
                                     pad=plotsize*5,
                                     fontsize=plotsize*5)
    if g == 2:
        for j in range(0, nc, step):
            axs[-1, j].set_xlabel(f'{logess[j]:2.0f}',
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
            ax.set_title(titles[c+1], pad=plotsize*9, fontsize=plotsize*5)
        if g == 2:
            ax.set_xticklabels(xticklabels, fontsize=ticklabels)
        axd.append(ax)
    axds.append(axd)

for g, folder in enumerate(folders):

    given = dfss[g][0].Given[0]
    lows = pd.pivot_table(dfss[g][0].loc[df.Time == ts[-1]],
                 values='a2Seenmean',
                 index=[rowindex],
                 columns=['logES']).sort_index(axis=0,
                                            ascending=False)
    lows = lows.to_numpy()
    T = my.fitness(highs, lows, given, AA, RR)
    R = my.fitness(highs, highs, given, AA, RR)
    P = my.fitness(lows, lows, given, AA, RR)
    S = my.fitness(lows, highs, given, AA, RR)

    Z = my.gamecolors(T, R, P, S)
    greys = np.full([*Z.shape], [0.8, 0.8, 0.8, 1.0])
    mask = (Z == my.colormap['white'])
    Z[mask] = greys[mask]

    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            y = [T[i, j], R[i, j], P[i, j], S[i, j]]
            for line in axss[g][i][j].get_lines():
                line.remove()
            axss[g][i][j].plot(xaxis,
                           y,
                           c=Z[i, j],
                           linewidth=3,
                           marker='o',
                           markerfacecolor='white',
                           markersize=plotsize/3)

    for j, trait in enumerate(traits):
        Z = pd.pivot_table(dfss[g].loc[df.Time == ts[-1]],
                           values=trait,
                           index=[rowindex],
                           columns=['logES']).sort_index(axis=0,
                                                ascending=False)
        axds[g][j].imshow(Z, vmin=0, vmax=traitvmaxs[j])

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
