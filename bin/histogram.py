#! /usr/bin/env python

from glob import glob
import imageio.v2 as iio
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

traits = ['a2Seen0', 'a2Seen31', 'a2Seen63']
folders = ['none', 'r']
movie = False

fslarge = 32 # Label font size
fssmall = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

dfs = []
for folder in folders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.frq'))), ignore_index=True)
    dfs.append(df)

df = dfs[0]
ts = df.Time.unique()
if movie:
    frames = []
else:
    ts = [ts[-1]]
given = df.Given[0]
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
xmin = logess[0]
xmax = logess[-1]
xlabel = 'Substitutability of $\it{B}$'
ymin = alphas[-1]
ymax = alphas[0]
ylabel = 'Value of $\it{B}$'
xticklabels = [round(xmin),
                round((xmin + xmax)/2),
                round(xmax)]
yticklabels = [round(ymin, 1),
                round((ymin + ymax)/2, 1),
                round(ymax, 1)]
everyx = int(nc/2)
everyy = int(nr/2)

fig = plt.figure(figsize=(16, 8))
fig.supxlabel(xlabel,
                x=0.525,
                y=0.03,
                fontsize=fslarge)
fig.supylabel(ylabel,
                x=0.05,
                y=0.52,
                fontsize=fslarge)

bins = [0.0, 0.5, 1.0]
my_cmap = plt.get_cmap('viridis')
colors = my_cmap(bins)

outergrid = fig.add_gridspec(1, 2, left=0.15, right=0.9, top=0.86, bottom=0.176)
innergrid = [outergrid[0].subgridspec(nrows=nr, ncols=nc, wspace=0, hspace=0), 
                outergrid[1].subgridspec(nrows=nr, ncols=nc, wspace=0, hspace=0)]
axss = [innergrid[0].subplots(), innergrid[1].subplots()]

for axs in axss:
    for axrow in axs:
        for ax in axrow:
            ax.set(xticks=[], yticks=[], xlim=(0.0, 1.0), ylim=(0.0, 1.0))
    for ax, loges in zip(axs[-1, ::everyx], logess[::everyx]):
        ax.set_xlabel(round(loges), fontsize=fssmall)
for ax, alpha in zip(axss[0][::everyy, 0], alphas[::everyy]):
    ax.set_ylabel(f'{alpha:1.1f}',
                    rotation='horizontal',
                    horizontalalignment='right',
                    verticalalignment='center',
                    fontsize=fssmall)

for t in ts:
    for df, axs in zip(dfs, axss):
        for axrow, alpha in zip(axs, alphas):
            for ax, loges in zip(axrow, logess):
                bottom = 0.0
                for trait, color in zip(traits, colors):
                    barheight = df.loc[(df['Time'] == t) & (df['logES'] > loges - 0.1) & (df['logES'] < loges + 0.1) & (df['alpha'] == alpha), trait].values[0]
                    ax.bar(x=0.5, height=barheight, bottom=bottom, width=1.0, color=color)
                    bottom = bottom + barheight
    if movie:
        text = fig.text(0.90,
                        0.90,
                        f't\n{t}',
                        fontsize=fssmall+4,
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
