#! /usr/bin/env python

from glob import glob
from math import log
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time

start_time = time.perf_counter ()

traits = ['ChooseGrainmean', 'MimicGrainmean']
traitlabels = ['Game types', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner']
folders = ['p', 'r']
alphafolders = ['alpha25', 'alpha50', 'alpha75']
alphas = [0.25, 0.50, 0.75]

letters = [['a', 'b', 'c'],
            ['d', 'e', 'f'],
            ['g', 'h', 'i']]

R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

def fitness(x, y, alpha):
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

dfss = []
for alphafolder in alphafolders:
    dfs = []
    for folder in folders:
        df = pd.concat(map(pd.read_csv, glob(os.path.join(alphafolder, 'discrete', folder, '*.csv'))), ignore_index=True)
        dfs.append(df)
    dfss.append(dfs)

df = dfss[0][0]
ts = df.Time.unique()
t = ts[-1]
givens = np.sort(pd.unique(df.Given))[::-1]
ess = np.sort(pd.unique(df.ES))

# Theory

b = a2max/a1max
rhos = 1.0 - 1.0/ess
nr = len(givens)
nc = len(rhos)
RR, GG = np.meshgrid(rhos, givens)
a20 = np.full([nc, nr], 0.0)
a21 = a20 + a2max/2.0
a22 = a20 + a2max

Zs = []

for alpha in alphas:
    T = fitness(a21, a20, alpha)
    R = fitness(a21, a21, alpha)
    P = fitness(a20, a20, alpha)
    S = fitness(a20, a21, alpha)
    Tu = fitness(a22, a21, alpha)
    Ru = fitness(a22, a22, alpha)
    Pu = fitness(a21, a21, alpha)
    Su = fitness(a21, a22, alpha)
    mask = (R < P)
    H = T[mask]
    T[mask] = S[mask]
    S[mask] = H
    H = R[mask]
    R[mask] = P[mask]
    P[mask] = H
    mask = (Su > R)
    T[mask] = Tu[mask]
    R[mask] = Ru[mask]
    P[mask] = Pu[mask]
    S[mask] = Su[mask]
    Z = np.full([nc, nr], 0.0)
    Z[(T < R) & (P < S)] = 0.9
    Z[(T >= R) & (P <= S)] = 0.5
    Z[(T > R) & (P > S)] = 0.1
    Zs.append(Z)

# Figure 

fslabel=24 # Label font size
fstick=16 # Tick font size

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig, axs = plt.subplots(nrows=len(alphas), ncols=len(traits)+1, figsize=(6*len(alphas), 6*(len(traits)+1)))

fig.supxlabel('Substitutability of $\it{A}$', x=0.52, y=0.03, fontsize=fslabel*1.35)
fig.supylabel('Partner\'s share of $\it{A}$', x=0.06, y=0.53, fontsize=fslabel*1.35, ha='center')

# Plots

minx = round(log(ess[0], 2))
maxx = round(log(ess[-1], 2))
miny = givens[-1]
maxy = givens[0]
xticklabels = [minx, round((minx + maxx)/2), maxx]
yticklabels = [miny, (miny + maxy)/2, maxy]

for axrow, letterrow in zip(axs, letters):
    for ax, letter, traitlabel in zip(axrow, letterrow, traitlabels):
        ax.text(0, nr*1.035, letter, fontsize=fslabel, weight='bold')
        ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
        if ax.get_subplotspec().is_first_row():
            ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 

extent = 0, nr, 0, nc

for axrow, Z, dfs in zip(axs, Zs, dfss):
    axrow[0].imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=1.0)
    for ax, df, trait in zip(axrow[1:], dfs, traits):
        df = df.loc[df.Time == t].copy()
        df[trait] = 1.0 - df[trait]
        df_piv = pd.pivot_table(df, values=trait, index=['Given'], columns=['ES']).sort_index(axis=0, ascending=False)
        ax.imshow(df_piv, extent=extent, cmap='magma', vmin=0, vmax=1.0)

plt.savefig('alphadiscrete.png', transparent=False)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
