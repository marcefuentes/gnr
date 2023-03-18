#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib import cm
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
titles = ['Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner']
folders = ['given100', 'given95', 'given50']
subfolders = ['none', 'p', 'r']

plotsize = 8

# Get data

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
t = ts[-1]
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1.0 - 1.0/pow(2.0, logess)
RR, AA = np.meshgrid(rhos, alphas)

# Figure properties

xlim=[0, 5]
ylim=[0.0, 2.0]
step = int(nr/2)
xaxis = [1, 2, 3, 4]
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

fig = plt.figure(figsize=(width*1.25, height))
fig.supxlabel(xlabel,
              x=0.502,
              y=0.04,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.04,
              y=0.5,
              fontsize=biglabels)
outergrid = fig.add_gridspec(nrows=len(folders),
                             ncols=len(traits),
                             left=0.22,
                             right=0.80)
axs = np.empty((len(folders),
                len(traits),
                nr,
                nc),
                dtype=object)

letter = ord('a')
letterposition = 4.8
for g, folder in enumerate(folders):
    for c, trait in enumerate(traits):
        grid = outergrid[g, c].subgridspec(nrows=nr,
                                           ncols=nc,
                                           wspace=0,
                                           hspace=0)
        axs[g, c] = grid.subplots()
        axs[g, c, 0, 0].set_title(chr(letter),
                                  fontsize=plotsize*5,
                                  pad = 10,
                                  weight='bold',
                                  loc='left')
        letter += 1
        for i, alpha in enumerate(alphas):
            for j, loges in enumerate(logess):
                axs[g, c, i, j].set(xticks=[], yticks=[])
                axs[g, c, i, j].set(xlim=xlim, ylim=ylim)
                for axis in ['top','bottom','left','right']:
                    axs[g, c, i, j].spines[axis].set_linewidth(0.1)
        if g == 0:
            axs[g, c, 0, 10].set_title(titles[c],
                                       pad=plotsize*9,
                                       fontsize=plotsize*5)
        for i in range(0, nr, step):
            axs[g, c, i, 0].set(yticks=[ylim[1]/2], yticklabels=[])
            if c == 0:
                axs[g, c, i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                                           rotation='horizontal',
                                           horizontalalignment='right',
                                           verticalalignment='center',
                                           y=0.3,
                                           fontsize=ticklabels)
        for j in range(0, nc, step):
            axs[g, c, -1, j].set(xticks=[xlim[1]/2], xticklabels=[])
            if folder == folders[-1]:
                axs[g, c, -1, j].set_xlabel(f'{logess[j]:2.0f}',
                                            x=0.3,
                                            fontsize=ticklabels)

df = dfsocial
m = df.Time == t
df = df.loc[m]
a2social = pd.pivot_table(df,
                          values='a2Seenmean',
                          index=['alpha'],
                          columns=['logES'])
a2social = a2social.sort_index(axis=0, ascending=False)
a2social = a2social.to_numpy()

for g, folder in enumerate(folders):
    for c, trait in enumerate(traits):

        df = dfs[g, 0]
        m = df.Time == t
        df = df.loc[m]
        given = df.Given.iloc[0]
        a2private = pd.pivot_table(df,
                                   values='a2Seenmean',
                                   index=['alpha'],
                                   columns=['logES'])
        a2private = a2private.sort_index(axis=0, ascending=False)
        a2private = a2private.to_numpy()

        T = my.fitness(a2social, a2private, given, AA, RR)
        R = my.fitness(a2social, a2social, given, AA, RR)
        P = my.fitness(a2private, a2private, given, AA, RR)
        S = my.fitness(a2private, a2social, given, AA, RR)

        df = dfs[g, c+1]
        m = df.Time == ts[-1]
        df = df.loc[m]
        Z = pd.pivot_table(df,
                           values=trait,
                           index=['alpha'],
                           columns=['logES'])
        Z = Z.sort_index(axis=0, ascending=False)
        Z = Z.to_numpy()

        for i, alpha in enumerate(alphas):
            for j, rho in enumerate(rhos):
                y = [T[i, j], R[i, j], P[i, j], S[i, j]]
                for line in axs[g, c, i, j].get_lines():
                    line.remove()
                axs[g, c, i, j].plot(xaxis,
                                     y,
                                     c=cm.viridis(Z[i, j]),
                                     linewidth=1,
                                     marker='o',
                                     markerfacecolor='white',
                                     markersize=plotsize/3)

# Save figure

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
