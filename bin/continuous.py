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

trait0s = ['nodilemmaRS', 'snowdrift', 'prisoner', 'prisonerRS']
traitlabel0s = ['Harmony game\n2$\it{R}$ < $\it{T}$ + $\it{S}$', 
                'Snowdrift', 
                'Prisoner\'s dilemma\n2$\it{R}$ > $\it{T}$ + $\it{S}$', 
                'Prisoner\'s dilemma\n2$\it{R}$ < $\it{T}$ + $\it{S}$']
traits = ['a2Seenmean', 'ChooseGrainmean', 'MimicGrainmean', 'wmean']
traitlabels = ['Effort to get $\it{B}$', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner', 'Fitness']
folders = ['none', 'p', 'r', 'pr', 'p8r', 'given0']

movie = False

filename = 'output'

fslabel = 32 # Label font size
fstick = 24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

dfgam = pd.concat(map(pd.read_csv, glob(os.path.join(folders[0], '*.gam'))), ignore_index=True)

dfs = []
for folder in folders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.csv'))), ignore_index=True)
    df.ChooseGrainmean = 1.0 - df.ChooseGrainmean
    df.MimicGrainmean = 1.0 - df.MimicGrainmean
    df['help'] = df.a2Seenmean*mymodule.R2*df.Given
    dfs.append(df)

df = dfs[0]
ts = df.Time.unique()
if movie:
    frames = []
else:
    ts = [ts[-1]]
logess = pd.unique(df.logES)
givens = pd.unique(df.Given)
alphas = pd.unique(df.alpha)

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
traitvmaxs = [mymodule.a2max, mymodule.a2max, mymodule.a2max, mymodule.fitness(np.array([mymodule.a2max]), np.array([mymodule.a2max]), np.array([0.0]), np.array([0.9]), np.array([5.0]))]
xticklabels = [round(xmin), round((xmin + xmax)/2), round(xmax)]
yticklabels = [round(ymin, 1), round((ymin + ymax)/2, 1), round(ymax, 1)]
extent = 0, nc, 0, nr

fig, axs = plt.subplots(nrows=len(folders)+1, ncols=len(traits), figsize=(6*len(traits), 6*(len(folders)+1)))
fig.supxlabel(xlabel, x=0.513, y=0.06, fontsize=fslabel*1.50)
fig.supylabel(ylabel, x=0.05, y=0.493, fontsize=fslabel*1.50, ha='center')

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
for ax, traitlabel in zip(axs[0], traitlabel0s):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)
for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

for t in ts:
    df = dfgam.loc[df.Time == t].copy()
    for ax, trait0 in zip(axs[0], trait0s):
        Z0 = pd.pivot_table(df, values=trait0, index=[rowindex], columns=['logES']).sort_index(axis=0, ascending=False)
        Z0 = Z0.to_numpy()
        Z = np.full([nr, nc, 4], mymodule.colormap[trait0])
        for i in range(Z0.shape[0]):
            for j in range(Z0.shape[1]):
                Z[i, j] = (1.0 - Z[i, j])*(1.0 - Z0[i, j]) + Z[i, j]
        ax.imshow(Z, extent=extent)
    for axrow, df in zip(axs[1:], dfs):
        for ax, trait, traitvmax in zip(axrow, traits, traitvmaxs):
            Z = pd.pivot_table(df.loc[df.Time == t], values=trait, index=[rowindex], columns=['logES']).sort_index(axis=0, ascending=False)
            ax.imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)
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
