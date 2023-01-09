#! /usr/bin/env python

from glob import glob
import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import pandas as pd
import time

start_time = time.perf_counter ()

traits0 = ['nodilemmaRS', 'snowdrift', 'prisoner', 'prisonerRS']
traitlabels0 = ['$\it{T}$ < $\it{R}$, $\it{P}$ < $\it{S}$\n2$\it{R}$ < $\it{T}$ + $\it{S}$', 
                'Snowdrift', 
                'Prisoner\'s dilemma', 
                'Prisoner\'s dilemma\n2$\it{R}$ < $\it{T}$ + $\it{S}$']
traitroots = ['a2Seen', 'ChooseGrain', 'MimicGrain', 'w']
traitlabels = ['Effort to get $\it{B}$', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner', 'Fitness']
folders = ['none', 'p', 'r', 'pr', 'p8r', 'given0']

movie = False

nbins = 64
filename = 'output'

fslabel = 32 # Label font size
fstick = 24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

dfgam = pd.concat(map(pd.read_csv, glob(os.path.join('none', '*.gam'))), ignore_index=True)

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
alphas = np.sort(pd.unique(df.alpha))[::-1]
nc = len(logess)
xmin = np.amin(logess)
xmax = np.amax(logess)
xlabel = 'Substitutability of $\it{B}$'
if len(givens) > 1:
    rows = givens
    ylabel = 'Partner\'s share of $\it{B}$'
    rowindex = 'Given'
if len(alphas) > 1:
    rows = alphas
    ylabel = 'Value of $\it{B}$'
    rowindex = 'alpha'
nr = len(rows)
ymin = np.amin(rows)
ymax = np.amax(rows)
xticklabels = [round(xmin), round((xmin + xmax)/2), round(xmax)]
yticklabels = [round(ymin, 1), round((ymin + ymax)/2, 1), round(ymax, 1)]
everyx = int(nc/2)
everyy = int(nr/2)
extent = 0, nc, 0, nr

fig = plt.figure(figsize=(len(traitroots)*6, (len(folders)+1)*6))
fig.supxlabel(xlabel, x=0.513, y=0.06, fontsize=fslabel*1.50)
fig.supylabel(ylabel, x=0.05, y=0.493, fontsize=fslabel*1.50, ha='center')

bins = np.arange(start=0, stop=nbins, step=1)
my_cmap = plt.get_cmap('magma')
colors = my_cmap(bins/nbins)

outergrid = fig.add_gridspec(len(folders)+1, len(traitroots), left=0.15, right=0.9, top=0.86, bottom=0.176)

axs = outergrid.subplots()
letter = ord('a')
for axrow in axs:
    for ax in axrow:
        if letter <= ord('z'): 
            ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
        else:
            ax.text(0, nr*1.035, 'a' + chr(letter - 26), fontsize=fslabel, weight='bold')
        ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
        if ax.get_subplotspec().is_first_row():
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
        letter += 1
for ax, traitlabel in zip(axs[0], traitlabels0):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)
for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

for t in ts:
    df = dfgam.loc[df.Time == t].copy()
    for ax, trait in zip(axs[0], traits0):
        Z0 = pd.pivot_table(df, values=trait, index=[rowindex], columns=['logES']).sort_index(axis=0, ascending=False)
        Z0 = Z0.to_numpy()
        Z = np.full([nr, nc, 4], mymodule.colormap[trait])
        for i in range(Z0.shape[0]):
            for j in range(Z0.shape[1]):
                Z[i, j] = (1.0 - Z[i, j])*(1.0 - Z0[i, j]) + Z[i, j]
        ax.imshow(Z, extent=extent)
    for outerrow, df in enumerate(dfs):
        for outercol, traitroot in enumerate(traitroots):
            traits = []
            for b in bins:
                traits.append(traitroot + str(b)) 
            innergrid = outergrid[outerrow+1, outercol].subgridspec(nrows=nr, ncols=nc, wspace=0, hspace=0)
            axs = innergrid.subplots()
            for axrow, row in zip(axs, rows):
                for ax, loges in zip (axrow, logess):
                    ax.set(xticks=[], yticks=[], xlim=(0.0, 1.0), ylim=(0.0, 1.0))
                    bottom = 0.0
                    for trait, color in zip(traits, colors):
                        barheight = df.loc[(df['Time'] == t) & (df['logES'] > loges - 0.1) & (df['logES'] < loges + 0.1) & (df[rowindex] == row), trait].values[0]
                        ax.bar(x=0.5, height=barheight, bottom=bottom, width=1.0, color=color)
                        bottom = bottom + barheight
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
