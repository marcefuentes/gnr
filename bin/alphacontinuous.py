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
alphafolders = ['alpha75', 'alpha50', 'alpha25']
alphas = [0.75, 0.50, 0.25]
datafolder = 'continuous/pr'
movie = False

numa2 = 32
filename = 'continuous'
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

fslabel = 32 # Label font size
fstick = 18 # Tick font size
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

dfs = []
for alphafolder in alphafolders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(alphafolder, datafolder, '*.csv'))), ignore_index=True)
    dfs.append(df)

df = dfs[0]
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
nc = len(rhos)
nr = len(givens)
minx = logess[0]
maxx = logess[-1]
miny = givens[-1]
maxy = givens[0]
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Partner\'s share of $\it{B}$'
pivindex = 'Given'

RRR, GGG = np.meshgrid(np.repeat(rhos, numa2), np.repeat(givens, numa2))

traitvmaxs = [a2max, a2max]
xticklabels = [round(minx), round((minx + maxx)/2), round(maxx)]
yticklabels = [round(miny, 1), round((miny + maxy)/2, 1), round(maxy, 1)]
extent = 0, nr, 0, nc
extenta2 = 0, nr*numa2, 0, nc*numa2
prisoner = [0.5, 0.0, 0.0, 1.0]
snowdrift = [0.0, 1.0, 1.0, 1.0]
nodilemma = [1.0, 1.0, 1.0, 1.0]
green = [0.0, 1.0, 0.0, 1.0]

b = a2max/a1max
X, Y = np.meshgrid(np.linspace(0.0, a2max, num=numa2), np.linspace(a2max, 0.0, num=numa2))
X = np.tile(A=X, reps=[nr, nc])
Y = np.tile(A=Y, reps=[nr, nc])

fig, axs = plt.subplots(nrows=len(alphafolders), ncols=len(traitlabels), figsize=(6*len(traitlabels), 6*len(alphafolders)))
fig.supxlabel(xlabel, x=0.512, y=0.03, fontsize=fslabel*1.25)
fig.supylabel(ylabel, x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

letter = ord('a')
for axrow in axs:
    for ax, traitlabel in zip(axrow, traitlabels):
        if ax.get_subplotspec().is_first_row():
            ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)
        if ax.get_subplotspec().is_first_col():
            ax.text(0, nr*numa2*1.035, chr(letter), fontsize=fslabel, weight='bold')
            ax.set(xticks=[0, numa2*nc/2, numa2*nc], yticks=[0, numa2*nr/2, numa2*nr], xticklabels=[]) 
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
        else:
            ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
            ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
        letter += 1

for axrow, alpha in zip(axs, alphas):
    AAA = np.full([nc*numa2, nr*numa2], alpha)
    T = fitness(Y, X, GGG, AAA, RRR)
    R = fitness(Y, Y, GGG, AAA, RRR)
    P = fitness(X, X, GGG, AAA, RRR)
    S = fitness(X, Y, GGG, AAA, RRR)
    mask = (R < P)
    H = T[mask]
    T[mask] = S[mask]
    S[mask] = H
    H = R[mask]
    R[mask] = P[mask]
    P[mask] = H
    Z = np.full((numa2*nr, numa2*nc, 4), green)
    Z[(T > R) & (P > S)] = prisoner
    Z[(T >= R) & (P <= S) & (R != P)] = snowdrift
    Z[((T < R) & (P < S)) | (R == P)] = nodilemma

    axrow[0].imshow(Z, extent=extenta2)

for t in ts:

    if movie:
        text = fig.text(0.90, 0.90, f't\n{t}', fontsize=fstick+4, color='grey', ha='right')

    for axrow, df in zip(axs, dfs):
        for ax, trait, traitvmax in zip(axrow[1:], traits, traitvmaxs):
            df = df.loc[df.Time == t].copy()
            df[trait] = 1.0 - df[trait]
            df_piv = pd.pivot_table(df, values=trait, index=[pivindex], columns=['logES']).sort_index(axis=0, ascending=False)
            ax.imshow(df_piv, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)

    if movie:
        plt.savefig('temp.png', transparent=False)
        frames.append(iio.imread('temp.png'))
        os.remove('temp.png')
        text.remove()
    else:
        plt.savefig(filename + '.png', transparent=False)

plt.close()

if movie:
    iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
