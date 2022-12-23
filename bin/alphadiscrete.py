#! /usr/bin/env python

from glob import glob
import os
import imageio.v2 as iio
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

movie = False

filename = 'alphadiscrete'
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

fslabel=36 # Label font size
fstick=24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def fitness(x, y, given, alpha, rho):
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - given) + x*R2*given
    w = q1*q2
    mask = (w > 0.0) & (rho == 0.0)
    w[mask] = pow(q1[mask], 1.0 - alpha[mask])*pow(q2[mask], alpha[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = (1.0 - alpha[mask])*pow(q1[mask], rho[mask]) + alpha[mask]*pow(q2[mask], rho[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = pow(w[mask], 1.0/rho[mask])
    mask = (rho > 0.0)
    w[mask] = pow((1.0 - alpha[mask])*pow(q1[mask], rho[mask]) + alpha[mask]*pow(q2[mask], rho[mask]), 1.0/rho[mask])
    return w

dfss = []
for alphafolder in alphafolders:
    dfs = []
    for folder in folders:
        df = pd.concat(map(pd.read_csv, glob(os.path.join(alphafolder, 'discrete', folder, '*.csv'))), ignore_index=True)
        dfs.append(df)
    dfss.append(dfs)

df = dfss[0][0]
ts = df.Time.unique()
if movie:
    frames = []
else:
    ts = [ts[-1]]
givens = np.sort(pd.unique(df.Given))[::-1]
if givens[0] > 0.9999999:
    givens[0] = 0.9999999
logess = np.sort(pd.unique(df.logES))
rhos = 1.0 - 1.0/pow(2.0, logess)
nr = len(givens)
nc = len(rhos)
minx = logess[0]
maxx = logess[-1]
xlabel = 'Substitutability of $\it{B}$'
miny = givens[-1]
maxy = givens[0]
ylabel = 'Partner\'s share of $\it{B}$'
pivindex = 'Given'

b = a2max/a1max
RR, GG = np.meshgrid(rhos, givens)

a20 = np.full([nc, nr], 0.0)
a21 = a20 + a2max/2.0
a22 = a20 + a2max

black = [0.2, 0.0, 0.2, 1.0]
cyan = [0.0, 1.0, 1.0, 1.0]
white = [1.0, 1.0, 1.0, 1.0]
green = [0.0, 1.0, 0.0, 1.0]

Zs = []

for alpha in alphas:
    T = fitness(a21, a20)
    R = fitness(a21, a21)
    P = fitness(a20, a20)
    S = fitness(a20, a21)
    Tu = fitness(a22, a21)
    Ru = fitness(a22, a22)
    Pu = fitness(a21, a21)
    Su = fitness(a21, a22)
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
    Z = np.full([nc, nr, 4], [0.0, 1.0, 0.0, 1.0])
    Z[(T < R) & (P < S)] = white
    Z[(T >= R) & (P <= S)] = cyan
    Z[(T > R) & (P > S)] = black
    Zs.append(Z)

fig, axs = plt.subplots(nrows=len(alphas), ncols=len(traits)+1, figsize=(6*len(alphas), 6*(len(traits)+1)))

fig.supxlabel('Substitutability of $\it{B}$', x=0.52, y=0.03, fontsize=fslabel*1.35)
fig.supylabel('Partner\'s share of $\it{B}$', x=0.06, y=0.53, fontsize=fslabel*1.35, ha='center')

extent = 0, nr, 0, nc

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

for axrow, Z in zip(axs, Zs):
    axrow[0].imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=1.0)
for t in ts:
    for axrow, dfs in zip(axs, dfss):
        for ax, df, trait in zip(axrow[1:], dfs, traits):
            df = df.loc[df.Time == t].copy()
            df[trait] = 1.0 - df[trait]
            df_piv = pd.pivot_table(df, values=trait, index=['Given'], columns=['logES']).sort_index(axis=0, ascending=False)
            ax.imshow(df_piv, extent=extent, cmap='magma', vmin=0, vmax=1.0)
    if movie:
        plt.savefig('temp.png', transparent=False)
        frames.append(iio.imread('temp.png'))
        os.remove('temp.png')
    else:
        plt.savefig(filename + '.png', transparent=False)

plt.close()

if movie:
    iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
