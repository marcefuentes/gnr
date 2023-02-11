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

traits = ['a2Seenmean',
            'ChooseGrainmean',
            'MimicGrainmean',
            'wmean']
traitlabels = ['Effort to get $\it{B}$',
                'Sensitivity for\nchoosing partner',
                'Sensitivity for\nmimicking partner',
                'Fitness']
folders = ['given0', 'none', 'p', 'r', 'pr', 'p8r']

movie = False

filename = 'output'

fslabel = 32 # Label font size
fstick = 24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

dfs = []
for folder in folders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.csv'))),
                    ignore_index=True)
    df.ChooseGrainmean = 1.0 - df.ChooseGrainmean
    df.MimicGrainmean = 1.0 - df.MimicGrainmean
    #df['help'] = df.a2Seenmean*mymodule.R2*df.Given
    dfs.append(df)

df = dfs[1]
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
    RR, GG = np.meshgrid(rhos, givens)
    ymin = givens[-1]
    ymax = givens[0]
    ylabel = 'Partner\'s share of $\it{B}$'
    rowindex = 'Given'
else:
    nr = len(alphas)
    GG = np.full([nr, nc], givens[0])
if len(alphas) > 1:
    nr = len(alphas)
    RR, AA = np.meshgrid(rhos, alphas)
    ymin = alphas[-1]
    ymax = alphas[0]
    ylabel = 'Value of $\it{B}$'
    rowindex = 'alpha'
else:
    nr = len(givens)
    AA = np.full([nr, nc], alphas[0])

traitvmaxs = [mymodule.a2max,
                mymodule.a2max,
                mymodule.a2max,
                mymodule.fitness(np.array([mymodule.a2max]),
                    np.array([mymodule.a2max]),
                    np.array([0.0]),
                    np.array([0.9]),
                    np.array([5.0]))]
xticklabels = [round(xmin),
                round((xmin + xmax)/2),
                round(xmax)]
yticklabels = [round(ymin, 1),
                round((ymin + ymax)/2, 1),
                round(ymax, 1)]
extent = 0, nc, 0, nr

zeros = np.zeros([nr, nc])
low = np.full([nr, nc], mymodule.a2low)
high = np.full([nr, nc], mymodule.a2high)
R = mymodule.fitness(high, high, GG, AA, RR)
P = mymodule.fitness(low, low, GG, AA, RR)

fig, axs = plt.subplots(nrows=len(folders)+1,
                        ncols=len(traits),
                        figsize=(6*len(traits), 6*(len(folders)+1)))
fig.supxlabel(xlabel,
                x=0.513,
                y=0.06,
                fontsize=fslabel*1.5)
fig.supylabel(ylabel,
                x=0.05,
                y=0.493,
                fontsize=fslabel*1.5,
                ha='center')

letter = ord('a')
for axrow in axs:
    for ax in axrow:
        if letter <= ord('z'): 
            ax.text(0,
                    nr*1.035,
                    chr(letter),
                    fontsize=fslabel,
                    weight='bold')
        else:
            ax.text(0,
                    nr*1.035,
                    'a' + chr(letter - 26),
                    fontsize=fslabel,
                    weight='bold')
        letter += 1
        ax.set(xticks=[0, nc/2, nc],
                yticks=[0, nr/2, nr],
                xticklabels=[],
                yticklabels=[])
        if ax.get_subplotspec().is_first_row():
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

a2eq = np.copy(zeros)
weq = np.copy(zeros)
xeq = np.copy(zeros)

T = mymodule.fitness(high, low, GG, AA, RR)
S = mymodule.fitness(low, high, GG, AA, RR)
Z = np.full([nr, nc, 4], mymodule.colormap['default'])
mymodule.gametypes(T, R, P, S, Z)
ax = axs[0, 0]
ax.imshow(Z, extent=extent)
ax.set_title('Game types', pad=50.0, fontsize=fslabel)

Z = np.zeros([nr, nc])
ax = axs[0, 1]
ax.imshow(Z, extent=extent, cmap='viridis', vmin=0, vmax=1.0)
ax.set_title('Prisoner\'s dilemma\n$\it{T}$ + $\it{S}$ - 2$\it{R}$',
                pad=50.0,
                fontsize=fslabel)

for t in ts:
    for axrow, df in zip(axs[1:], dfs):
        for ax, trait, traitvmax in zip(axrow, traits, traitvmaxs):
            Z = pd.pivot_table(df.loc[df.Time == t],
                                values=trait,
                                index=[rowindex],
                                columns=['logES']).sort_index(axis=0,
                                                            ascending=False)
            ax.imshow(Z,
                    extent=extent,
                    cmap='viridis',
                    vmin=0,
                    vmax=traitvmax)
    if movie:
        text = fig.text(0.90,
                        0.93,
                        f't\n{t}',
                        fontsize=fstick+4,
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
