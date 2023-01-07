#! /usr/bin/env python

from glob import glob
import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time

start_time = time.perf_counter ()

traits = ['a2Seen0', 'a2Seen31', 'a2Seen63']
folders = ['none', 'p']
movie = False

filename = 'histogram'
fslabel = 32 # Label font size
fstick = 18 # Tick font size
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
givens = np.sort(pd.unique(df.Given))[::-1]
logess = np.sort(pd.unique(df.logES))
rhos = 1.0 - 1.0/pow(2.0, logess)
alphas = np.sort(pd.unique(df.alpha))[::-1]
nc = len(rhos)
xmin = logess[0]
xmax = logess[-1]
xlabel = 'Substitutability of $\it{B}$'

if len(givens) > 1:
    nr = len(givens)
    ymin = givens[-1]
    ymax = givens[0]
    ylabel = 'Partner\'s share of $\it{B}$'
    rowname = 'Given'
    rows = givens
else:
    nr = len(alphas)
if len(alphas) > 1:
    nr = len(alphas)
    ymin = alphas[-1]
    ymax = alphas[0]
    ylabel = 'Value of $\it{B}$'
    rowname = 'alpha'
    rows = alphas
else:
    nr = len(givens)

everyx = int(nc/2)
everyy = int(nr/2)
xticklabels = [round(xmin), round((xmin + xmax)/2), round(xmax)]
yticklabels = [round(ymin, 1), round((ymin + ymax)/2, 1), round(ymax, 1)]

fig = plt.figure(figsize=(16, 8))
fig.supxlabel(xlabel, x=0.525, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.05, y=0.52, fontsize=fslabel)

outergrid = fig.add_gridspec(1, 2, left=0.15, right=0.9, top=0.86, bottom=0.176)

for t in ts:
    for df, outer in zip(dfs, outergrid):
        innergrid = outer.subgridspec(nrows=nr, ncols=nc, wspace=0, hspace=0)
        axs = innergrid.subplots()
        for axrow in axs:
            for ax in axrow:
                ax.set(xticks=[], yticks=[], ylim=(0.0, 1.1))
        for ax, loges in zip(axs[-1, ::everyx], logess[::everyx]):
            ax.set_xlabel(round(loges), fontsize=fstick)
        for ax, row in zip(axs[::everyy, 0], rows[::everyy]):
            ax.set_ylabel(f'{row:1.1f}', rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)
        for axrow, row in zip(axs, rows):
            for ax, loges in zip(axrow, logess):
                for b in ['a2Seen0', 'a2Seen31', 'a2Seen63']:
                    barheight = df.loc[(df['Time'] == t) & (df['logES'] == loges) & (df[rowname] == row), b]
                    ax.bar(x=b, height=barheight, align='edge')

    if movie:
        text = fig.text(0.90, 0.90, f't\n{t}', fontsize=fstick+4, color='grey', ha='right')
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
