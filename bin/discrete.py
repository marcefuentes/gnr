#! /usr/bin/env python

from glob import glob
from math import log
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time

start_time = time.perf_counter ()

letters = [['a', 'b', 'c', 'd', 'e'],
            ['f', 'g', 'h', 'i', 'j'],
            ['k', 'l', 'm', 'n', 'o'],
            ['p', 'q', 'r', 's', 't'],
            ['u', 'v', 'w', 'x', 'y'],
            ['z', 'aa', 'ab', 'ac', 'ad']]

#traits = ['a2Seensd', 'help', 'wsd', 'ChooseGrainsd', 'MimicGrainsd']
traits = ['a2Seenmean', 'help', 'wmean', 'ChooseGrainmean', 'MimicGrainmean']
traitlabels = ['Effort to get $\it{A}$', 'Help', 'Fitness', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner']
traitvmaxs = [0.5, 1.0, 1.0, 0.6, 0.6]
folders = ['none', 'p', 'r', 'pr', 'p8r']

alpha = 0.5
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

def fitness(x, y):
    q1 = (1.0 - y)*R1
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

dfs = []
for folder in folders:
    dfs.append(pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.csv'))), ignore_index=True))

t = dfs[0].Time.iat[-1]

dfts = []
for df in dfs:
    dfts.append(df.loc[df.Time == t].copy())
for df in dfts:
    df['help'] = df[traits[0]]*R2*df.Given

R = R2/R1
b = a2max/a1max
ess = np.sort(pd.unique(dfts[0].ES))
rhos = 1.0 - 1.0/ess
givens = np.sort(pd.unique(dfts[0].Given))[::-1]
nc = len(rhos)
nr = len(givens)
givens[0] = 0.9999999
RR, GG = np.meshgrid(rhos, givens)
a2Ds = np.full([nc, nr], 0.0)
a2Cs = np.full([nc, nr], a2max/2.0)
Ts = fitness(a2Cs, a2Ds)
Rs = fitness(a2Cs, a2Cs)
Ps = fitness(a2Ds, a2Ds)
Ss = fitness(a2Ds, a2Cs)
xeqss = np.full([nc, nr], 0.0)
mask = (Ts < Rs) & (Ps < Ss)
xeqss[mask] = 1.0
mask = (Ts > Rs) & (Ps < Ss) & (Rs - Ss - Ts + Ps != 0.0)
xeqss[mask] = (Ps[mask] - Ss[mask])/(Rs[mask] - Ss[mask] - Ts[mask] + Ps[mask])
weqss = (Ts + Ss)*xeqss*(1.0 - xeqss) + Rs*xeqss*xeqss + Ps*(1.0 - xeqss)*(1.0 - xeqss)
a2eqss = xeqss*a2max/2.0
helpeqss = a2eqss*R2*GG 
Zs = [a2eqss, helpeqss, weqss, np.ones([nc, nr])*0.06, np.ones([nc, nr])*0.06]

fslabel=36 # Label font size
fstick=24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig, axs = plt.subplots(nrows=len(folders)+1, ncols=len(traits), figsize=(6*len(traits), 6*(len(folders)+1)))

fig.supxlabel('Substitutability of $\it{A}$', x=0.513, y=0.05, fontsize=fslabel*1.25)
fig.supylabel('Partner\'s share of $\it{A}$', x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

extent = 0, nr, 0, nc

# Top row of subplots

for ax, Z, traitvmax in zip(axs[0], Zs, traitvmaxs):
    ax.imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)

# Remaining rows of plots

for axrow, df in zip(axs[1:], dfts):
    for ax, trait, traitvmax in zip(axrow, traits, traitvmaxs):
        if 'Grain' in trait: df[trait] = 1.0 - df[trait]
        df_piv = pd.pivot_table(df, values=trait, index=['Given'], columns=['ES']).sort_index(axis=0, ascending=False)
        ax.imshow(df_piv, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)

# All plots

for axrow in axs:
    for ax in axrow:
        ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])

givens[0] = 1.0
minx = round(log(ess[0], 2))
maxx = round(log(ess[-1], 2))
miny = givens[-1]
maxy = givens[0]

xticklabels = [minx, round((minx + maxx)/2), maxx]
yticklabels = [miny, (miny + maxy)/2, maxy]
for axrow in axs:
    axrow[0].set_yticklabels(yticklabels, fontsize=fstick) 
for ax in axs[-1]:
    ax.set_xticklabels(xticklabels, fontsize=fstick)

for ax, traitlabel in zip(axs[0], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

for axrow, letterrow in zip(axs, letters):
    for ax, letter in zip(axrow, letterrow):
        ax.text(0, nr*1.035, letter, fontsize=fslabel, weight='bold')

plt.savefig('discrete.png', transparent=False)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')