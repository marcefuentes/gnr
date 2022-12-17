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

traits = ['ChooseGrainmean', 'MimicGrainmean']
traitlabels = ['Game types', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner']
folders = ['p', 'r']
alphafolders = ['alpha25', 'alpha50', 'alpha75']
alphas = [0.25, 0.50, 0.75]

movie = False
if movie:
    frames = []

npoints = 128
filename = 'alphacontinuous'
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

# Figure

fslabel=24 # Label font size
fstick=16 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

letters = [['a', 'b', 'c'],
            ['d', 'e', 'f'],
            ['g', 'h', 'i']]

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

# Simulations

dfss = []
for alphafolder in alphafolders:
    dfs = []
    for folder in folders:
        df = pd.concat(map(pd.read_csv, glob(os.path.join(alphafolder, 'continuous', folder, '*.csv'))), ignore_index=True)
        dfs.append(df)
    dfss.append(dfs)

df = dfss[0][0]
ts = df.Time.unique()
if movie == False: ts = [ts[-1]]
givens = np.sort(pd.unique(df.Given))[::-1]
ess = np.sort(pd.unique(df.ES))
givens[0] = 0.9999999
rhos = 1.0 - 1.0/ess
nr = len(givens)
nc = len(rhos)
mingiven = givens[-1]
maxgiven = givens[0]

# Theory

b = a2max/a1max
Rq = R2/R1
a2x = np.linspace(0.0, a2max, num=npoints)
a2y = np.linspace(a2max, 0.0, num=npoints)
X, Y = np.meshgrid(a2x, a2y)
RR, GG = np.meshgrid(rhos, givens)
TT = b*Rq*(1.0 - GG)

minx = round(log_ess[0])
maxx = round(log_ess[-1])
miny = mingiven
maxy = maxgiven

xticklabels = [minx, round((minx + maxx)/2), maxx]
yticklabels = [miny, (miny + maxy)/2, maxy]
extent = 0, num, 0, num
extentZ = 0, npoints*num, 0, npoints*num

Zsss = []
for alpha in alphas:
    QQ = Rq*pow(TT*alpha/(1.0 - alpha), 1.0/(RR - 1.0))
    Z = np.zeros([npoints, npoints])

    A = np.full([npoints, npoints], alpha)
    Zss = np.empty((0, npoints*num))
    for given in givens:
        G = np.full([npoints, npoints], given)
        Zs = np.empty((npoints, 0))
        for rho in rhos:
            Z = np.zeros([npoints, npoints])
            Rh = np.full([npoints, npoints], rho)
            T = fitness(Y, X, G, A, Rh)
            R = fitness(Y, Y, G, A, Rh)
            P = fitness(X, X, G, A, Rh)
            S = fitness(X, Y, G, A, Rh)
            mask = (R < P)
            H = R[mask]
            R[mask] = P[mask]
            P[mask] = H
            H = T[mask]
            T[mask] = S[mask]
            S[mask] = H
            Z[(T > R) & (P > S)] = 0.1
            Z[(T >= R) & (P <= S)] = 0.5
            Z[(T < R) & (P < S)] = 0.9
            Z[(T == R) & (P == R)] = 1.0
            #Z = np.tril(Z, k=-1)
            Z = np.ma.masked_where(Z == 0.0, Z)
            Zs = np.append(Zs, Z, axis=1)
        Zss = np.append(Zss, Zs, axis=0)
    Zsss.append(Zss)

fig, axs = plt.subplots(nrows=len(folders)+1, ncols=len(traits), figsize=(6*len(traits), 6*(len(folders)+1)))
fig.supxlabel('Substitutability of $\it{B}$', x=0.513, y=0.05, fontsize=fslabel*1.25)
fig.supylabel('Partner\'s share of $\it{B}$', x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

# Plots

minx = round(log(ess[0], 2))
maxx = round(log(ess[-1], 2))
xticklabels = [minx, round((minx + maxx)/2), maxx]
yticklabels = [0.0, 0.5, 1.0]
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

for axrow, Zss in zip(axs, Zsss):
    axrow[0].imshow(Zss, extent=extent, cmap='magma', vmin=0, vmax=1.0)
for t in ts:
    for axrow, dfs in zip(axs, dfss):
        for ax, df, trait in zip(axrow[1:], dfs, traits):
            df = df.loc[df.Time == t].copy()
            df[trait] = 1.0 - df[trait]
            df_piv = pd.pivot_table(df, values=trait, index=['Given'], columns=['ES']).sort_index(axis=0, ascending=False)
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
