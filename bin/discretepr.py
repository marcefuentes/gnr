#! /usr/bin/env python

from glob import glob
import imageio.v2 as iio
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
folders = ['high', 'mid', 'low']
subfolder = 'pr'

movie = False
rows = folders
plotsize = 4

dfs = []
for folder in folders:
    filelist = glob(os.path.join(folder, subfolder, '*.csv'))
    df = pd.concat(map(pd.read_csv, filelist),
                   ignore_index=True)
    for trait in traits:
        df[trait] = 1.0 - df[trait]
    dfs.append(df)

df = dfs[0]
ts = df.Time.unique()
if movie:
    frames = []
else:
    ts = [ts[-1]]
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

fig, axs = plt.subplots(nrows=len(rows),
                        ncols=len(titles),
                        figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.513,
              y=0.03,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.04,
              y=0.493,
              fontsize=biglabels)

for ax in fig.get_axes():
    ax.set(xticks=xticks, yticks=yticks)
    ax.set(xticklabels=[], yticklabels=[])
    ax.text(0,
            letterposition,
            chr(letter),
            transform=ax.transAxes,
            fontsize=plotsize*5,
            weight='bold')
    letter += 1
for i, row in enumerate(rows):
    axs[i, 0].set_yticklabels(yticklabels, fontsize=ticklabels)
for j, title in enumerate(titles):
    axs[0, j].set_title(title, pad=plotsize*9, fontsize=plotsize*5)
    axs[-1, j].set_xticklabels(xticklabels,
                               x=0.47,
                               fontsize=ticklabels)

for t in ts:

    for g, folder in enumerate(folders):

        df = dfss[g][0]
        m = [df.Time == t
        df = df.loc[m]
        given = df.Given.iloc[0]
        lows = pd.pivot_table(df,
                              values='a2low',
                              index=[rowindex],
                              columns=['logES'])
        lows = lows.sort_index(axis=0, ascending=False)
        lows = lows.to_numpy()
        highs = pd.pivot_table(df,
                               values='a2high',
                               index=[rowindex],
                               columns=['logES'])
        highs = highs.sort_index(axis=0, ascending=False)
        highs = highs.to_numpy()
        T = my.fitness(highs, lows, given, AA, RR)
        R = my.fitness(highs, highs, given, AA, RR)
        P = my.fitness(lows, lows, given, AA, RR)
        S = my.fitness(lows, highs, given, AA, RR)

        Z = my.gamecolors(T, R, P, S)
        axs[g, 0].imshow(Z)

        for j, trait in enumerate(traits):
            df = dfss[g][j]
            m = df.Time == t
            df = df.loc[m]
            Z = pd.pivot_table(df,
                               values=trait,
                               index=[rowindex],
                               columns=['logES'])
            Z = Z.sort_index(axis=0, ascending=False)
            axs[g, j + 1].imshow(Z, vmin=0, vmax=traitvmaxs[j])

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