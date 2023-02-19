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

trait0labels = ['Games (lower)',
                '$\it{R}$ - $\it{P}$',
                'Games (upper)',
                '$\it{T}$ + $\it{S}$ - 2$\it{R}$']
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

num = 512

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
    dfs.append(df)

df = dfs[1]
ts = df.Time.unique()
if movie:
    frames = []
else:
    ts = [ts[-1]]
given = df.Given[0]
if given > 0.9999999:
    given = 0.9999999
alphas = np.sort(pd.unique(df.alpha))[::-1]
rowindex = 'alpha'
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
alphas = np.linspace(alphas[0], alphas[-1], num=num)
logess = np.linspace(logess[0], logess[-1], num=num)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
MRT0 = mymodule.b*mymodule.Rq
Q0 = mymodule.Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)
MRT = MRT0*(1.0 - given)
Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
a2lows = [0.0, a2eq]
a2highs = [(1.0 + 4.0*a2social)/5.0, mymodule.a2max]

xmin = logess[0]
xmax = logess[-1]
xlabel = 'Substitutability of $\it{B}$'
ymin = alphas[-1]
ymax = alphas[0]
ylabel = 'Value of $\it{B}$'

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
extent2 = 0, num, 0, num

fig, axs = plt.subplots(nrows=len(folders)+1,
                        ncols=len(traits),
                        figsize=(6*len(traits), 6*(len(folders)+1)))
fig.supxlabel(xlabel, x=0.513, y=0.06, fontsize=fslabel*1.5)
fig.supylabel(ylabel, x=0.05, y=0.493, fontsize=fslabel*1.5)

letter = ord('a')
for axrow in axs:
    for ax in axrow:
        if letter <= ord('z'): 
            textl = chr(letter)
        else:
            textl = 'a' + chr(letter - 26)
        letter += 1
        if ax.get_subplotspec().is_first_row():
            pixels = num
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        else:
            pixels = nc
        ax.set(xticks=[0, pixels/2, pixels],
                yticks=[0, pixels/2, pixels],
                xticklabels=[],
                yticklabels=[])
        ax.text(0, pixels*1.035, textl, fontsize=fslabel, weight='bold')
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
for ax, traitlabel in zip(axs[0], trait0labels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)
for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

for i, (a2low, a2high) in enumerate(zip(a2lows, a2highs)):

    low = np.full([num, num], a2low)
    high = np.full([num, num], a2high)
    T = mymodule.fitness(high, low, given, AA, RR)
    R = mymodule.fitness(high, high, given, AA, RR)
    P = mymodule.fitness(low, low, given, AA, RR)
    S = mymodule.fitness(low, high, given, AA, RR)
    Z = np.full([num, num, 4], mymodule.colormap['red'])
    mymodule.gamecolors(T, R, P, S, Z)
    axs[0, 2*i].imshow(Z, extent=extent2)

    if i == 0:
        Z = np.zeros([num, num])
        mask = (R > P)
        Z[mask] = R[mask] - P[mask]
        Z = np.ma.masked_where(Z == 0.0, Z)
        cmap = plt.cm.viridis
        cmap.set_bad(color='white')
        axs[0, 2*i+1].imshow(Z, extent=extent2, cmap=cmap)
    else:
        Z = np.zeros([num, num])
        mask = mymodule.prisoner(T, R, P, S) | (mymodule.deadlock(T, R, P, S) & (2.0*P < T + S)) | (R == P)
        Z[mask] = 1.0 - (2.0*R[mask] - T[mask] - S[mask])
        Z = np.ma.masked_where(Z == 0.0, Z)
        cmap = plt.cm.viridis
        cmap.set_bad(color='white')
        axs[0, 2*i+1].imshow(Z, extent=extent2, cmap=cmap)

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
                        0.90,
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
