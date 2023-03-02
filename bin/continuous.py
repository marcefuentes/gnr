#! /usr/bin/env python

from glob import glob
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import pandas as pd
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

titlegs = ['Games',
           '$\it{R}$ - $\it{P}$',
           '$\it{T}$ + $\it{S}$ - 2$\it{R}$']
traits = ['a2Seenmean',
          'ChooseGrainmean',
          'MimicGrainmean',
          'wmean']
titles = ['Effort to get $\it{B}$',
               'Sensitivity for\nchoosing partner',
               'Sensitivity for\nmimicking partner',
               'Fitness']
traitvmaxs = [mymodule.a2max,
              mymodule.a2max,
              mymodule.a2max,
              mymodule.fitness(np.array([mymodule.a2max]),
                               np.array([mymodule.a2max]),
                               np.array([0.0]),
                               np.array([0.9]),
                               np.array([5.0]))]
folders = ['given0', 'none', 'p', 'r', 'pr', 'p8r']

movie = False
rows = folders + 1
plotsize = 4

dfs = []
for folder in folders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.csv'))),
                    ignore_index=True)
    df.ChooseGrainmean = 1.0 - df.ChooseGrainmean
    df.MimicGrainmean = 1.0 - df.MimicGrainmean
    dfs.append(df)

df = dfs[1]
ts = df.Time.unique()
if movie:
    frames = []
else:
    ts = [ts[-1]]
given = df.Given[0]
alphas = np.sort(pd.unique(df.alpha))[::-1]
rowindex = 'alpha'
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1.0 - 1.0/pow(2.0, logess)
RR, AA = np.meshgrid(rhos, alphas)

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
              x=0.515,
              y=0.06,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.03,
              y=0.493,
              fontsize=biglabels)

outergrid = fig.add_gridspec(nrows=1),
                             ncols=len(titlegs),
                             left=0.15,
                             right=0.85,
                             top=0.8,
                             bottom=0.2)

axsg = []
for g, title in enumerate(titlegs):
    grid = outergrid[0, g].subgridspec(nrows=nr,
                                       ncols=nc,
                                       wspace=0,
                                       hspace=0)
    axs = grid.subplots()
    axs[0, int(num/2)].set_title(title,
                                 pad=plotsize*5,
                                 fontsize=plotsize*5)
    axs[0, 0].set_title(chr(letter),
                        fontsize=plotsize*5,
                        weight='bold',
                        loc='left')
    letter += 1

    for ax in fig.get_axes():
        ax.set(xticks=[], yticks=[])
        ax.set(xticklabels=[])
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(0.1)
    if g == 0:
        for i in range(0, num, step):
            axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                                 rotation='horizontal',
                                 horizontalalignment='right',
                                 verticalalignment='center',
                                 fontsize=ticklabels)
    axsg.append(axs)

#grid = fig.add_gridspec(nrows=len(folders),
#                        ncols=len(titlegs),
#                        left=0.15,
#                        right=0.85,
#                        top=0.8,
#                        bottom=0.2)
#
#axs = grid.subplots()

for g, title in enumerate(titles):
    
    for ax in fig.get_axes():
        ax.set(xticks=xticks, yticks=yticks)
        ax.set(xticklabels=[], yticklabels=[])
        if letter <= ord('z'): 
            text = chr(letter)
        else:
            text = 'a' + chr(letter - 26)
        ax.text(0,
                letterposition,
                text,
                transform=ax.transAxes,
                fontsize=plotsize*5,
                weight='bold')
        letter += 1
    for i, row in enumerate(rows):
        axs[i, 0].set_yticklabels(yticklabels, fontsize=ticklabels)
    for j, title in enumerate(titlegs):
        ax = axs[0, j]
        ax.set_title(title, pad=plotsize*10, fontsize=plotsize*5)
        pos = ax.get_position()
        newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
        ax.set_position(newpos)
        
    for j, title in enumerate(titles):
        axs[1, j].set_title(title, pad=plotsize*10, fontsize=plotsize*5)
        axs[-1, j].set_xticklabels(xticklabels, fontsize=ticklabels)

for i, title in enumerate(titlegs):
    for i, alpha in enumerate(alphas):
        AA = np.full([ext, ext], alpha)
        for j, rho in enumerate(rhos):

            xmin = 0.0
            xmax = mymodule.a2eq(given, alpha, rho)
            ymin = mymodule.a2eq(0.0, alpha, rho)
            ymax = mymodule.a2max
            x = np.linspace(xmin, xmax, num=ext)
            y = np.linspace(ymax, ymin, num=ext)
            X, Y = np.meshgrid(x, y)
            RR = np.full([ext, ext], rho)
            T = mymodule.fitness(Y, X, given, AA, RR)
            R = mymodule.fitness(Y, Y, given, AA, RR)
            P = mymodule.fitness(X, X, given, AA, RR)
            S = mymodule.fitness(X, Y, given, AA, RR)

            Z = mymodule.gamecolors(T, R, P, S)
            axss[0][i][j].imshow(Z, extent=extent)

            Z = R - P
            axss[1][i][j].imshow(Z, extent=extent, vmin=-1, vmax=1)

            Z = T + S - 2.0*R
            mask = R < P
            Z[mask] = T[mask] + S[mask] - 2.0*P[mask]
            axss[2][i][j].imshow(Z, extent=extent, vmin=-1, vmax=1)

for t in ts:
    for i, df in enumerate(dfs):
        for j, trait in enumerate(traits):
            Z = pd.pivot_table(df.loc[df.Time == t],
                               values=trait,
                               index=[rowindex],
                               columns=['logES']).sort_index(axis=0,
                                                    ascending=False)
            axs[i, j].imshow(Z, vmin=0, vmax=traitvmaxs[j])
    if movie:
        text = fig.text(0.90,
                        0.93,
                        f't\n{t}',
                        fontsize=biglabels,
                        color='grey',
                        ha='right')
        plt.savefig('temp.png', transparent=False)
        text.remove()
        frames.append(iio.imread('temp.png'))
        os.remove('temp.png')
    else:
        plt.savefig(filename + '.png', transparent=False)

plt.close()

if movie:
    iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
