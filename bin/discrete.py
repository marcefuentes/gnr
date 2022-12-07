#! /usr/bin/env python

from glob import glob
from math import log
import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time

start_time = time.perf_counter ()

traits = ['a2Seenmean', 'help', 'wmean', 'ChooseGrainmean', 'MimicGrainmean']
traitlabels = ['Effort to get $\it{A}$', 'Help', 'Fitness', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner']
traitvmaxs = [1.0, 2.0, 1.5, 1.0, 1.0]
folders = ['none', 'p', 'r', 'pr', 'p8r']

letters = [['a', 'b', 'c', 'd', 'e'],
            ['f', 'g', 'h', 'i', 'j'],
            ['k', 'l', 'm', 'n', 'o'],
            ['p', 'q', 'r', 's', 't'],
            ['u', 'v', 'w', 'x', 'y'],
            ['z', 'aa', 'ab', 'ac', 'ad'],
            ['ae', 'af', 'ag', 'ah', 'ai']]

movie = False

alpha = 0.5
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

def fitness(x, y):
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - GG) + x*R2*GG
    w = q1*q2
    mask = (w > 0.0) & (RR == 0.0)
    w[mask] = pow(q1[mask], alpha)*pow(q2[mask], 1.0 - alpha)
    mask = (w > 0.0) & (RR < 0.0)
    w[mask] = alpha*pow(q1[mask], RR[mask]) + (1.0 - alpha)*pow(q2[mask], RR[mask])
    mask = (w > 0.0) & (RR < 0.0)
    w[mask] = pow(w[mask], 1.0/RR[mask])
    mask = (RR > 0.0)
    w[mask] = pow(alpha*pow(q1[mask], RR[mask]) + (1.0 - alpha)*pow(q2[mask], RR[mask]), 1.0/RR[mask])
    return w

# Simulations

dfs = []
for folder in folders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.csv'))), ignore_index=True)
    df.ChooseGrainmean = 1.0 - df.ChooseGrainmean
    df.MimicGrainmean = 1.0 - df.MimicGrainmean
    df['help'] = df.a2Seenmean*R2*df.Given
    dfs.append(df)

df = dfs[0]
ts = df.Time.unique()
if movie == False: ts = [ts[-1]]
givens = np.sort(pd.unique(df.Given))[::-1]
ess = np.sort(pd.unique(df.ES))

# Theory

b = a2max/a1max
rhos = 1.0 - 1.0/ess
nr = len(givens)
nc = len(rhos)
RR, GG = np.meshgrid(rhos, givens)
a20 = np.full([nc, nr], 0.0)
a2 = a20
w = fitness(a20, a20)

for c in range(2):
    a20 = a20 + a2max*c/2.0
    a21 = a20 + a2max/2.0
    T = fitness(a21, a20)
    R = fitness(a21, a21)
    P = fitness(a20, a20)
    S = fitness(a20, a21)
    mask = (T < R) & (P < S)
    a2[mask] = a2max*(1.0 + c)/2.0
    w[mask] = R[mask]
    mask = (T >= R) & (P < S) & (R - S - T + P != 0.0)
    T = T[mask]
    R = R[mask]
    P = P[mask]
    S = S[mask]
    x = (P - S)/(R - S - T + P)
    a2[mask] = x*a2max*(1.0 + c)/2.0
    w[mask] = (T + S)*x*(1.0 - x) + R*x*x + P*(1.0 - x)*(1.0 - x)

helps = a2*R2*GG 
Zs = [a2, helps, w, a20, a20]

# Figure

fslabel=36 # Label font size
fstick=24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
frames = []

fig, axs = plt.subplots(nrows=len(folders)+1, ncols=len(traits), figsize=(6*len(traits), 6*(len(folders)+1)))

fig.supxlabel('Substitutability of $\it{A}$', x=0.513, y=0.05, fontsize=fslabel*1.25)
fig.supylabel('Partner\'s share of $\it{A}$', x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

extent = 0, nr, 0, nc
givens[0] = 1.0
minx = round(log(ess[0], 2))
maxx = round(log(ess[-1], 2))
miny = givens[-1]
maxy = givens[0]
xticklabels = [minx, round((minx + maxx)/2), maxx]
yticklabels = [miny, (miny + maxy)/2, maxy]

for axrow in axs:
    for ax in axrow:
        ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
    axrow[0].set_yticklabels(yticklabels, fontsize=fstick) 
for ax in axs[-1]:
    ax.set_xticklabels(xticklabels, fontsize=fstick)
for ax, traitlabel in zip(axs[0], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)
for axrow, letterrow in zip(axs, letters):
    for ax, letter in zip(axrow, letterrow):
        ax.text(0, nr*1.035, letter, fontsize=fslabel, weight='bold')

for t in ts:
    # Row 0: theory (no cooperation)
    for ax, Z, traitvmax in zip(axs[0], Zs, traitvmaxs):
        ax.imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)
    # Remaining rows: simulations 
    for axrow, df in zip(axs[1:], dfs):
        for ax, trait, traitvmax in zip(axrow, traits, traitvmaxs):
            df_piv = pd.pivot_table(df.loc[df.Time == t], values=trait, index=['Given'], columns=['ES']).sort_index(axis=0, ascending=False)
            ax.imshow(df_piv, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)
    plt.savefig('output.png', transparent=False)
    if movie:
        frames.append(iio.imread('output.png'))
        os.remove('output.png')
if movie:
    iio.mimsave('output.gif', frames)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
