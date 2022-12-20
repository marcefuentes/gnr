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

traits = ['a2Seenmean', 'ChooseGrainmean', 'MimicGrainmean', 'wmean']
traitlabels = ['Effort to get $\it{B}$', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner', 'Fitness']
traitvmaxs = [1.0, 1.0, 1.0, 1.8]
folders = ['none', 'p', 'r', 'pr', 'p8r']
given = 0.95

movie = False

numa2 = 64
filename = 'output'
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

# Figure

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

dfs = []
for folder in folders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.csv'))), ignore_index=True)
    df.ChooseGrainmean = 1.0 - df.ChooseGrainmean
    df.MimicGrainmean = 1.0 - df.MimicGrainmean
    df['help'] = df.a2Seenmean*R2*given
    dfs.append(df)

df = dfs[0]
ts = df.Time.unique()
if movie:
    frames = []
else:
    ts = [ts[-1]]
alphas = np.sort(pd.unique(df.alpha))[::-1]
ess = np.sort(pd.unique(df.ES))
rhos = 1.0 - 1.0/pow(2.0, ess)
nr = len(alphas)
nc = len(rhos)

b = a2max/a1max
Rq = R2/R1
MRT0 = b*Rq
RR, AA = np.meshgrid(rhos, alphas)
X, Y = np.meshgrid(np.linspace(0.0, a2max, num=numa2), np.linspace(a2max, 0.0, num=numa2))
X = np.tile(A=X, reps=[nr, nc])
Y = np.tile(A=Y, reps=[nr, nc])
RRR, AAA = np.meshgrid(np.repeat(rhos, numa2), np.repeat(alphas, numa2))

minx = ess[0]
maxx = ess[-1]
miny = alphas[-1]
maxy = alphas[0]

xticklabels = [minx, round((minx + maxx)/2), maxx]
yticklabels = [miny, (miny + maxy)/2, maxy]
extent = 0, nr, 0, nc
extenta2 = 0, numa2*nr, 0, numa2*nc
prisoner = np.full((numa2*nr, numa2*nc, 4), [0.5, 0.0, 0.0, 1.0])
snowdrift = np.full((numa2*nr, numa2*nc, 4), [0.0, 1.0, 1.0, 1.0])
nodilemma = np.full((numa2*nr, numa2*nc, 4), [1.0, 1.0, 1.0, 1.0])

Z = np.full((numa2*nr, numa2*nc, 4), [0.0, 1.0, 0.0, 1.0])
T = fitness(Y, X, given, AAA, RRR)
R = fitness(Y, Y, given, AAA, RRR)
P = fitness(X, X, given, AAA, RRR)
S = fitness(X, Y, given, AAA, RRR)
mask = (R < P)
H = R[mask]
R[mask] = P[mask]
P[mask] = H
H = T[mask]
T[mask] = S[mask]
S[mask] = H
mask = (T > R) & (P > S)
Z[mask] = prisoner[mask]
mask = (T >= R) & (P <= S) & (R != P)
Z[mask] = snowdrift[mask]
mask = ((T < R) & (P < S)) | (R == P)
Z[mask] = nodilemma[mask]
#Z = np.tril(Z, k=-1)
#Z = np.ma.masked_where(Z == 0.0, Z)

MRT = MRT0*(1.0 - given)
Q0 = Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2socialss = a2max/(1.0 + Q0*b)
Q = Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2eqss = a2max/(1.0 + Q*b)

Mss = [[a2socialss, a2socialss*R2*given, fitness(a2socialss, a2socialss, given, AA, RR)], [a2eqss, a2eqss*R2*given, fitness(a2eqss, a2eqss, given, AA, RR)]]

fig, axs = plt.subplots(nrows=len(folders)+1, ncols=len(traits), figsize=(6*len(traits), 6*(len(folders)+1)))
fig.delaxes(axs[0, 1])
fig.delaxes(axs[0, 2])
fig.delaxes(axs[0, 3])
fig.supxlabel('Value of $\it{B}$', x=0.513, y=0.05, fontsize=fslabel*1.25)
fig.supylabel('Partner\'s share of $\it{B}$', x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

letter = ord('b')
for axrow in axs:
    for ax in axrow:
        ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
        if ax.get_subplotspec().is_first_row():
            ax.set(xticks=[0, numa2*nr/2, numa2*nc], yticks=[0, numa2*nr/2, numa2*nc], xticklabels=[], yticklabels=yticklabels)
            ax.set_title('Game types', pad=50.0, fontsize=fslabel)
            ax.text(0, numa2*nr*1.035, 'a', fontsize=fslabel, weight='bold')
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        else:
            ax.set(xticks=[0, nc/2, nr], yticks=[0, nc/2, nr], xticklabels=[], yticklabels=[])
            ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
            letter += 1
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 

for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

axs[0, 0].imshow(Z, extent=extenta2)

for t in ts:

    if movie:
        text = fig.text(0.80, 0.80, f't\n{t}', fontsize=fstick+4, color='grey', ha='right')

    for axrow, df in zip(axs[1:], dfs):
        for ax, trait, traitvmax in zip(axrow, traits, traitvmaxs):
            df_piv = pd.pivot_table(df.loc[df.Time == t], values=trait, index=['alpha'], columns=['ES']).sort_index(axis=0, ascending=False)
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
