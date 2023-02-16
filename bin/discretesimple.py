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
                '$\it{T}$ + $\it{S}$ - 2*$\it{R}$',
                'Sensitivity for\nchoosing partner',
                'Sensitivity for\nmimicking partner']
folders = ['a2init75', 'a2init50', 'a2init25']
subfolder = 'pr'
a2lows = [0.50, 0.25, 0.00]
given = 0.95

num = 1001    # Number of subplot rows and columns
filename = 'output'

movie = False

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
nr = len(alphas)
nc = len(logess)
alphas = np.linspace(alphas[0], alphas[-1], num=num)
logess = np.linspace(logess[0], logess[-1], num=num)
rhos = 1.0 - 1.0/pow(2.0, logess)
RR, AA = np.meshgrid(rhos, alphas)

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
extent2 = 0, num, 0, num

fig, axs = plt.subplots(nrows=len(a2lows),
                        ncols=len(traitlabels),
                        figsize=(6*len(traitlabels), 6*(len(a2lows))))
fig.supxlabel(xlabel, x=0.513, y=0.02, fontsize=fslabel*1.3)
fig.supylabel(ylabel, x=0.03, y=0.493, fontsize=fslabel*1.3)

letter = ord('a')
for axrow in axs:
    for ax, traitlabel in zip(axrow, traitlabels):
        if traitlabel == 'Games':
            ax.text(0, 
                    num*1.035,
                    chr(letter),
                    fontsize=fslabel*0.8,
                    weight='bold')
        else:
            ax.text(0, 
                    nr*1.035,
                    chr(letter),
                    fontsize=fslabel*0.8,
                    weight='bold')
        letter += 1
        if ax.get_subplotspec().is_first_col():
            ax.set(xticks=[0, num/2, num],
                    yticks=[0, num/2, num],
                    xticklabels=[])
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
        else:
            ax.set(xticks=[0, nc/2, nc],
                            yticks=[0, nr/2, nr],
                            xticklabels=[],
                            yticklabels=[])
        if ax.get_subplotspec().is_first_row():
            ax.set_title(traitlabel, pad=40.0, fontsize=fslabel)
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)

for axrow, a2low in zip(axs, a2lows):

    low = np.full([num, num], a2low)
    high = np.full([num, num], a2low + 0.5)
    T = mymodule.fitness(high, low, given, AA, RR)
    R = mymodule.fitness(high, high, given, AA, RR)
    P = mymodule.fitness(low, low, given, AA, RR)
    S = mymodule.fitness(low, high, given, AA, RR)
    Z = np.full([num, num, 4], mymodule.colormap['red'])
    mymodule.gamecolors(T, R, P, S, Z)
    axrow[0].imshow(Z, extent=extent2)

    Z = np.zeros([num, num])
    mask = mymodule.prisoner(T, R, P, S) & (2.0*R > T + S)
    Z[mask] = 1.0 - (2.0*R[mask] - T[mask] - S[mask])
    Z = np.ma.masked_where(Z == 0.0, Z)
    cmap = plt.cm.viridis
    cmap.set_bad(color='white')
    axrow[1].imshow(Z, extent=extent, cmap=cmap, vmin=0, vmax=traitvmaxs[0])

for t in ts:
    for axrow, df in zip(axs, dfs):
        for ax, trait, traitvmax in zip(axrow[2:], traits, traitvmaxs):
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
