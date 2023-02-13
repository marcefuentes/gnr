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

traits = ['ChooseGrainmean',
            'MimicGrainmean']
traitlabels = ['Games',
                'Sensitivity for\nchoosing partner',
                'Sensitivity for\nmimicking partner']
folders = ['given95', 'given50']
folders = ['given95']
givens = [0.95, 0.5]
subfolder = 'pr'

movie = False

filename = 'output'

fslabel = 32 # Label font size
fstick = 24 # Tick font size
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

alphas = np.sort(pd.unique(df.alpha))[::-1]
rowindex = 'alpha'
logess = np.sort(pd.unique(df.logES))
rhos = 1.0 - 1.0/pow(2.0, logess)
nr = len(alphas)
nc = len(logess)
low = np.full([nr, nc], mymodule.a2low)
high = np.full([nr, nc], mymodule.a2high)
RR, AA = np.meshgrid(rhos, alphas)
R = mymodule.fitness(high, high, 0.0, AA, RR)
P = mymodule.fitness(low, low, 0.0, AA, RR)

xmin = logess[0]
xmax = logess[-1]
xlabel = 'Substitutability of $\it{B}$'
ymin = alphas[-1]
ymax = alphas[0]
ylabel = 'Value of $\it{B}$'

traitvmaxs = [mymodule.a2max,
                mymodule.a2max]
xticklabels = [round(xmin),
                round((xmin + xmax)/2),
                round(xmax)]
yticklabels = [round(ymin, 1),
                round((ymin + ymax)/2, 1),
                round(ymax, 1)]
extent = 0, nc, 0, nr

fig, axs = plt.subplots(nrows=len(givens),
                        ncols=len(traitlabels),
                        figsize=(6*len(traitlabels), 6*(len(givens))))
fig.supxlabel(xlabel, x=0.513, y=0.01, fontsize=fslabel*1.2)
fig.supylabel(ylabel, x=0.03, y=0.493, fontsize=fslabel*1.2)

letter = ord('a')
for axrow in axs:
    for ax, traitlabel in zip(axrow, traitlabels):
        ax.text(0,
                nr*1.035,
                chr(letter),
                fontsize=fslabel,
                weight='bold')
        ax.set(xticks=[0, nc/2, nc],
                yticks=[0, nr/2, nr],
                xticklabels=[],
                yticklabels=[])
        letter += 1
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
        if ax.get_subplotspec().is_first_row():
            ax.set_title(traitlabel, pad=40.0, fontsize=fslabel*0.9)
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)

for axrow, given in zip(axs, givens):

    T = mymodule.fitness(high, low, given, AA, RR)
    S = mymodule.fitness(low, high, given, AA, RR)
    Z = np.full([nr, nc, 4], mymodule.colormap['default'])
    mymodule.gametypes(T, R, P, S, Z)
    axrow[0].imshow(Z, extent=extent)

for t in ts:
    for axrow, df in zip(axs, dfs):
        for ax, trait, traitvmax in zip(axrow[1:], traits, traitvmaxs):
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
