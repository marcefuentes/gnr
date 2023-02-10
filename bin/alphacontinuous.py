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
traits = ['ChooseGrainmean', 'MimicGrainmean']
traitlabels = ['Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner']
alphafolders = ['alpha75', 'alpha50', 'alpha25']
alphas = [0.75, 0.50, 0.25]
datafolderg = 'none'
datafolders = ['p', 'r']
movie = False

filename = 'output'

fslabel = 32 # Label font size
fstick = 24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

dfgs = []
dfss = []
for alphafolder in alphafolders:
    dfg = pd.concat(map(pd.read_csv, glob(os.path.join(alphafolder, datafolderg, '*.gam'))), ignore_index=True)
    dfgs.append(dfg)
    dfs = []
    for datafolder in datafolders:
        dfs.append(pd.concat(map(pd.read_csv, glob(os.path.join(alphafolder, datafolder, '*.csv'))), ignore_index=True))
    dfss.append(dfs)

df = dfss[0][0]
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
xticklabels = [round(xmin), round((xmin + xmax)/2), round(xmax)]
yticklabels = [round(ymin, 1), round((ymin + ymax)/2, 1), round(ymax, 1)]
extent = 0, nc, 0, nr

fig, axs = plt.subplots(nrows=len(alphafolders), ncols=len(trait0s)+len(traits), figsize=(6*(len(trait0s)+len(traits)), 6*len(alphafolders)))
fig.supxlabel(xlabel, x=0.512, y=0.02, fontsize=fslabel*1.25)
fig.supylabel(ylabel, x=0.07, y=0.493, fontsize=fslabel*1.25, ha='center')

alltraitlabels = traitlabel0s + traitlabels
letter = ord('a')
for axrow in axs:
    for ax, alltraitlabel in zip(axrow, alltraitlabels):
        ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
        ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
        if ax.get_subplotspec().is_first_row():
            ax.set_title(alltraitlabel, pad=50.0, fontsize=fslabel)
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
        letter += 1

for t in ts:
    for axrow, dfg, dfs in zip(axs, dfgs, dfss):
        df = dfg.loc[dfg.Time == t].copy()
        for ax, trait0 in zip(axrow[:len(trait0s)], trait0s):
            Z0 = pd.pivot_table(df, values=trait0, index=[rowindex], columns=['logES']).sort_index(axis=0, ascending=False)
            Z0 = Z0.to_numpy()
            Z = np.full([nr, nc, 4], mymodule.colormap[trait0])
            for i in range(Z0.shape[0]):
                for j in range(Z0.shape[1]):
                    Z[i, j] = (1.0 - Z[i, j])*(1.0 - Z0[i, j]) + Z[i, j]
            ax.imshow(Z, extent=extent)
        for ax, df, trait in zip(axrow[len(trait0s):], dfs, traits): 
            df = df.loc[df.Time == t].copy()
            df[trait] = 1.0 - df[trait]
            Z = pd.pivot_table(df, values=trait, index=[rowindex], columns=['logES']).sort_index(axis=0, ascending=False)
            ax.imshow(Z, extent=extent, cmap='viridis', vmin=0, vmax=1)
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
