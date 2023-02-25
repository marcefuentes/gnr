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
folders = ['a2init75', 'a2init50', 'a2init25']
subfolder = 'pr'

movie = False

fslarge = 32 # Label font size
fssmall = 24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

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
letterposition = nr*1.035
xmin = logess[0]
xmax = logess[-1]
ymin = alphas[-1]
ymax = alphas[0]
xticks = [0, nc/2, nc]
yticks = [0, nr/2, nr]
xticklabels = [f'{xmin:2.0f}',
                f'{(xmin + xmax)/2.0:2.0f}',
                f'{xmax:2.0f}']
yticklabels = [f'{ymin:3.1f}',
                f'{(ymin + ymax)/2.0:3.1f}',
                f'{ymax:3.1f}']
extent= 0, nc, 0, nr
cmap = plt.cm.viridis
cmap.set_bad(color='white')
traitvmaxs = [mymodule.a2max,
                mymodule.a2max]

fig, axs = plt.subplots(nrows=len(folders),
                        ncols=len(titles),
                        figsize=(6*len(titles), 6*len(dfs)))
fig.supxlabel(xlabel,
                x=0.513,
                y=0.02,
                fontsize=fslarge*1.5)
fig.supylabel(ylabel,
                x=0.05,
                y=0.493,
                fontsize=fslarge*1.5)

for i, folder in enumerate(folders):
    for j, title in enumerate(titles):
        ax = axs[i, j]
        ax.set(xticks=xticks,
                yticks=yticks,
                xticklabels=[],
                yticklabels=[])
        ax.text(0,
                letterposition,
                chr(letter),
                fontsize=fslarge*0.8,
                weight='bold')
        letter += 1
        if ax.get_subplotspec().is_first_row():
            ax.set_title(title, pad=40.0, fontsize=fslarge)
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fssmall)
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fssmall)

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
    Z = np.full([nr, nc, 4], mymodule.colormap['red'])
    mymodule.gamecolors(T, R, P, S, Z)
    axs[i, 0].imshow(Z, extent=extent)

    Z = np.zeros([nr, nc])
    mask = mymodule.dilemma(T, R, P, S)
    Z[mask] = R[mask] - P[mask] + 0.000001
    Z = np.ma.masked_where(Z == 0.0, Z)
    axs[i, 1].imshow(Z, extent=extent, cmap=cmap, vmin=-1, vmax=1)

    Z = np.zeros([nr, nc])
    mask = mymodule.dilemma(T, R, P, S)
    Z[mask] = T[mask] + S[mask] - 2.0*R[mask]
    Z = np.ma.masked_where(Z == 0.0, Z)
    axs[i, 2].imshow(Z, extent=extent, cmap=cmap, vmin=-1, vmax=1)

for t in ts:
    for i, df in enumerate(dfs):
        for j, (trait, vmax) in enumerate(zip(traits, traitvmaxs)):
            Z = pd.pivot_table(df.loc[df.Time == t],
                        values=trait,
                        index=[rowindex],
                        columns=['logES']).sort_index(axis=0,
                                                    ascending=False)
            axs[i, j + 3].imshow(Z,
                                extent=extent,
                                vmin=0,
                                vmax=vmax)
    if movie:
        text = fig.text(0.90,
                        0.93,
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

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
