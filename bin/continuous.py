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

titles = ['Games (lower)',
                '$\it{R}$ - $\it{P}$',
                'Games (upper)',
                '$\it{T}$ + $\it{S}$ - 2$\it{R}$']
traits = ['a2Seenmean',
            'ChooseGrainmean',
            'MimicGrainmean',
            'wmean']
titletraits = ['Effort to get $\it{B}$',
                'Sensitivity for\nchoosing partner',
                'Sensitivity for\nmimicking partner',
                'Fitness']
folders = ['given0', 'none', 'p', 'r', 'pr', 'p8r']

movie = False

num = 1001

fslarge = 32 # Label font size
fssmall = 24 # Tick font size
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
rhos = 1.0 - 1.0/pow(2.0, logess)
RR, AA = np.meshgrid(rhos, alphas)
MRT0 = mymodule.b*mymodule.Rq
Q0 = mymodule.Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)
MRT = MRT0*(1.0 - given)
Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
a2lows = [0.01*a2eq, 0.99*a2eq]
a2highs = [0.99*a2social + 0.01*mymodule.a2max,
            0.01*a2social + 0.99*mymodule.a2max]

xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
letterposition = num*1.035
letterpositionr = nr*1.035
xmin = logess[0]
xmax = logess[-1]
ymin = alphas[-1]
ymax = alphas[0]
xticks = [0, num/2, num]
yticks = [0, num/2, num]
xticksr = [0, nc/2, nc]
yticksr = [0, nr/2, nr]
xticklabels = [f'{xmin:2.0f}',
                f'{(xmin + xmax)/2.0:2.0f}',
                f'{xmax:2.0f}']
yticklabels = [f'{ymin:3.1f}',
                f'{(ymin + ymax)/2.0:3.1f}',
                f'{ymax:3.1f}']
extentnum = 0, num, 0, num
extentr= 0, nc, 0, nr
cmap = plt.cm.viridis
cmap.set_bad(color='white')
traitvmaxs = [mymodule.a2max,
                mymodule.a2max,
                mymodule.a2max,
                mymodule.fitness(np.array([mymodule.a2max]),
                                    np.array([mymodule.a2max]),
                                    np.array([0.0]),
                                    np.array([0.9]),
                                    np.array([5.0]))]

fig, axs = plt.subplots(nrows=len(folders)+1,
                        ncols=len(traits),
                        figsize=(6*len(traits), 6*(len(folders)+1)))
fig.supxlabel(xlabel,
                x=0.513,
                y=0.06,
                fontsize=fslarge*1.5)
fig.supylabel(ylabel,
                x=0.05,
                y=0.493,
                fontsize=fslarge*1.5)

for i in range(axs):
    for j in range(traits):
        ax = axs[i, j]
        if letter <= ord('z'): 
            textl = chr(letter)
        else:
            textl = 'a' + chr(letter - 26)
        letter += 1
        if ax.get_subplotspec().is_first_row():
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
            ax.set(xticks=xticks,
                    yticks=yticks,
                    xticklabels=[],
                    yticklabels=[])
            ax.text(0,
                    letterposition,
                    textl,
                    fontsize=fslarge,
                    weight='bold')
        else:
            ax.set(xticks=xticksr,
                    yticks=yticksr,
                    xticklabels=[],
                    yticklabels=[])
            ax.text(0,
                    letterpositionr,
                    textl,
                    fontsize=fslarge,
                    weight='bold')
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fssmall) 
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fssmall)
for j in enumerate(titles):
    axs[0, j].set_title(title, pad=50.0, fontsize=fslarge)
for j in enumerate(titletraits):
    axs[1, j].set_title(title, pad=50.0, fontsize=fslarge)

for i, (a2low, a2high) in enumerate(zip(a2lows, a2highs)):

    low = np.full([num, num], a2low)
    high = np.full([num, num], a2high)
    T = mymodule.fitness(high, low, given, AA, RR)
    R = mymodule.fitness(high, high, given, AA, RR)
    P = mymodule.fitness(low, low, given, AA, RR)
    S = mymodule.fitness(low, high, given, AA, RR)
    Z = np.full([num, num, 4], mymodule.colormap['red'])
    mymodule.gamecolors(T, R, P, S, Z)
    axs[0, 2*i].imshow(Z, extent=extentnum)

    if i == 0:
        Z = np.zeros([num, num])
        mask = mymodule.dilemma(T, R, P, S)
        Z[mask] = R[mask] - P[mask]
        Z = np.ma.masked_where(Z == 0.0, Z)
        axs[0, 2*i+1].imshow(Z, extent=extentnum, cmap=cmap)
    else:
        Z = np.zeros([num, num])
        mask = mymodule.dilemma(T, R, P, S)
        Z[mask] = 1.0 - (2.0*R[mask] - T[mask] - S[mask])
        Z = np.ma.masked_where(Z == 0.0, Z)
        axs[0, 2*i+1].imshow(Z, extent=extentnum, cmap=cmap)

for t in ts:
    for axrow, df in zip(axs[1:], dfs):
        for ax, trait, traitvmax in zip(axrow, traits, traitvmaxs):
            Z = pd.pivot_table(df.loc[df.Time == t],
                                values=trait,
                                index=[rowindex],
                                columns=['logES']).sort_index(axis=0,
                                                            ascending=False)
            ax.imshow(Z,
                    extent=extentr,
                    vmin=0,
                    vmax=traitvmax)
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

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
