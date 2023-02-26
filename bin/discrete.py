#! /usr/bin/env python

from glob import glob
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import pandas as pd
import time

start_time = time.perf_counter ()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

titles = ['Games',
          '$\it{R}$ - $\it{P}$',
          '$\it{T}$ + $\it{S}$ - 2$\it{R}$',
          'Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner']
traits = ['ChooseGrainmean',
          'MimicGrainmean']
traitvmaxs = [mymodule.a2max,
              mymodule.a2max]
folders = ['a2init75', 'a2init50', 'a2init25']
subfolder = 'pr'

movie = False
rows = folders
plotsize = 4

dfs = []
for folder in folders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, subfolder, '*.csv'))),
                    ignore_index=True)
    df.ChooseGrainmean = 1.0 - df.ChooseGrainmean
    df.MimicGrainmean = 1.0 - df.MimicGrainmean
    dfs.append(df)

df = dfs[0]
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
xmin = logess[0]
xmax = logess[-1]
ymin = alphas[-1]
ymax = alphas[0]
xticks = [0, nc/2-0.5, nc-1]
yticks = [0, nr/2-0.5, nr-1]
xticklabels = [f'{xmin:2.0f}',
               f'{(xmin + xmax)/2.0:2.0f}',
               f'{xmax:2.0f}']
yticklabels = [f'{ymax:3.1f}',
               f'{(ymin + ymax)/2.0:3.1f}',
               f'{ymin:3.1f}']
width = plotsize*len(titles)
height = plotsize*len(rows)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig, axs = plt.subplots(nrows=len(rows),
                        ncols=len(titles),
                        figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.513,
              y=0.02,
              fontsize=plotsize*6+width/2)
fig.supylabel(ylabel,
              x=0.05,
              y=0.493,
              fontsize=plotsize*6+width/2)

for ax in fig.get_axes():
    ax.set(xticks=xticks,
           yticks=yticks,
           xticklabels=[],
           yticklabels=[])
    ax.text(0,
            letterposition,
            chr(letter),
            transform=ax.transAxes,
            fontsize=plotsize*5,
            weight='bold')
    letter += 1
for i, row in enumerate(rows):
    axs[i, 0].set_yticklabels(yticklabels, fontsize=plotsize*4)
for j, title in enumerate(titles):
    axs[0, j].set_title(title, pad=plotsize*10, fontsize=plotsize*6)
    axs[-1, j].set_xticklabels(xticklabels, fontsize=plotsize*4)

for i, df in enumerate(dfs):

    lows = pd.pivot_table(df.loc[df.Time == 1],
                values='a2low',
                index=[rowindex],
                columns=['logES']).sort_index(axis=0,
                                            ascending=False)
    lows = lows.to_numpy()
    highs = pd.pivot_table(df.loc[df.Time == 1],
                values='a2high',
                index=[rowindex],
                columns=['logES']).sort_index(axis=0,
                                            ascending=False)
    highs = highs.to_numpy()
    T = mymodule.fitness(highs, lows, given, AA, RR)
    R = mymodule.fitness(highs, highs, given, AA, RR)
    P = mymodule.fitness(lows, lows, given, AA, RR)
    S = mymodule.fitness(lows, highs, given, AA, RR)

    Z = mymodule.gamecolors(T, R, P, S)
    axs[i, 0].imshow(Z)

    N = mymodule.nodilemmacolors(T, R, P, S)

    Z = R - P
    axs[i, 1].imshow(Z, vmin=-1, vmax=1)
    axs[i, 1].imshow(N)

    Z = np.zeros([nr, nc])
    Z = T + S - 2.0*R
    axs[i, 2].imshow(Z, vmin=-1, vmax=1)
    axs[i, 2].imshow(N)

for t in ts:
    for i, df in enumerate(dfs):
        for j, trait in enumerate(traits):
            Z = pd.pivot_table(df.loc[df.Time == t],
                               values=trait,
                               index=[rowindex],
                               columns=['logES']).sort_index(axis=0,
                                                    ascending=False)
            axs[i, j + 3].imshow(Z, vmin=0, vmax=traitvmaxs[j])
    if movie:
        text = fig.text(0.90,
                        0.93,
                        f't\n{t}',
                        fontsize=width*2,
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

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
