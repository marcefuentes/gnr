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
folders = ['none', 'pr']

movie = False

filename = 'output'
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

def gametypes(a2c, a2d):
    mask = (mask0 & (T < R) & (P < S))
    Z[mask] = nodilemma[mask]
    a2eq[mask] = a2c[mask]
    weq[mask] = R[mask]
    mask = (mask0 & (T >= R) & (P <= S) & (R != P))
    Z[mask] = snowdrift[mask]
    xeq[mask] = (P[mask] - S[mask])/(R[mask] - S[mask] - T[mask] + P[mask])
    a2eq[mask] = a2c[mask]*xeq[mask] + a2d[mask]*(1.0 - xeq[mask])
    weq[mask] = (T[mask] + S[mask])*xeq[mask]*(1.0 - xeq[mask]) + R[mask]*xeq[mask]*xeq[mask] + P[mask]*(1.0 - xeq[mask])*(1.0 - xeq[mask])
    mask = (mask0 & (T > R) & (P > S))
    Z[mask] = prisoner[mask]
    a2eq[mask] = a2d[mask]
    weq[mask] = P[mask]
    pass

dfs = []
for folder in folders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.csv'))), ignore_index=True)
    df.ChooseGrainmean = 1.0 - df.ChooseGrainmean
    df.MimicGrainmean = 1.0 - df.MimicGrainmean
    df['help'] = df.a2Seenmean*R2*df.Given
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
ess = np.sort(pd.unique(df.ES))
rhos = 1.0 - 1.0/ess
alphas = np.sort(pd.unique(df.alpha))[::-1]
nc = len(rhos)
minx = int(log(ess[0], 2.0))
maxx = int(log(ess[-1], 2.0))
xlabel = 'Substitutability of $\it{B}$'

if len(givens) > 1:
    nr = len(givens)
    RR, GG = np.meshgrid(rhos, givens)
    miny = round(givens[-1], 1)
    maxy = round(givens[0], 1)
    ylabel = 'Partner\'s share of $\it{B}$'
    pivindex = 'Given'
else:
    nr = len(alphas)
    GG = np.full([nc, nr], givens[0])
if len(alphas) > 1:
    nr = len(alphas)
    RR, AA = np.meshgrid(rhos, alphas)
    miny = round(alpha[-1], 1)
    maxy = round(alpha[0], 1)
    ylabel = 'Value of $\it{B}$'
    pivindex = 'alpha'
else:
    nr = len(givens)
    AA = np.full([nc, nr], alphas[0])

xticklabels = [minx, round((minx + maxx)/2), maxx]
yticklabels = [miny, (miny + maxy)/2, maxy]
extent = 0, nr, 0, nc
prisoner = np.full((nr, nc, 4), [0.5, 0.0, 0.0, 1.0])
snowdrift = np.full((nr, nc, 4), [0.0, 1.0, 1.0, 1.0])
nodilemma = np.full((nr, nc, 4), [1.0, 1.0, 1.0, 1.0])
green = np.full((nr, nc, 4), [0.0, 1.0, 0.0, 1.0])

b = a2max/a1max
zeros = np.zeros([nr, nc])
a20 = np.copy(zeros)
a21 = np.full([nr, nc], a2max/2.0)
a22 = np.full([nr, nc], a2max)
w00 = fitness(a20, a20, zeros, AA, RR)
w11 = fitness(a21, a21, zeros, AA, RR)
w22 = fitness(a22, a22, zeros, AA, RR)
a2social = np.copy(zeros)
mask = (w11 > w00)
a2social[mask] = a21[mask]
mask = (w22 > w11)
a2social[mask] = a22[mask]
wsocial = fitness(a2social, a2social, zeros, AA, RR)

Z = np.copy(green)
a2eq = np.copy(zeros)
weq = np.copy(zeros)
xeq = np.copy(zeros)
w01 = fitness(a20, a21, GG, AA, RR)
w10 = fitness(a21, a20, GG, AA, RR)
w12 = fitness(a21, a22, GG, AA, RR)
w21 = fitness(a22, a21, GG, AA, RR)
w02 = fitness(a20, a22, GG, AA, RR)
w20 = fitness(a22, a20, GG, AA, RR)

mask0 = (wsocial == w00)
T = np.copy(w01)
R = np.copy(w00)
P = np.copy(w11)
S = np.copy(w10)
gametypes(a20, a21)

mask0 = (wsocial == w11) & (w20 > w22)
T = np.copy(w10)
R = np.copy(w11)
P = np.copy(w00)
S = np.copy(w01)
gametypes(a21, a20)

mask0 = (wsocial == w11) & (w20 <= w22)
T = np.copy(w12)
R = np.copy(w11)
P = np.copy(w22)
S = np.copy(w21)
gametypes(a21, a22)

mask0 = (wsocial == w22) & (w10 < w11)
T = np.copy(w21)
R = np.copy(w22)
P = np.copy(w11)
S = np.copy(w12)
gametypes(a22, a21)

mask0 = (wsocial == w22) & (w10 >= w11)
T = np.copy(w10)
R = np.copy(w11)
P = np.copy(w00)
S = np.copy(w01)
gametypes(a21, a20)

Mss = [[a2social, a2social*R2*GG, wsocial], [a2eq, a2eq*R2*GG, weq]]

fig, axs = plt.subplots(nrows=len(folders)+1, ncols=len(traits), figsize=(6*len(traits), 6*(len(folders)+1)))
fig.delaxes(axs[0, 1])
fig.delaxes(axs[0, 2])
fig.delaxes(axs[0, 3])
fig.supxlabel(xlabel, x=0.513, y=0.05, fontsize=fslabel*1.25)
fig.supylabel(ylabel, x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

letter = ord('b')
for axrow in axs:
    for ax in axrow:
        ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
        if ax.get_subplotspec().is_first_row():
            ax.set_title('Game types', pad=50.0, fontsize=fslabel)
            ax.text(0, nr*1.035, 'a', fontsize=fslabel, weight='bold')
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        else:
            ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
            letter += 1
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 

for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

axs[0, 0].imshow(Z, extent=extent)

for t in ts:

    if movie:
        text = fig.text(0.80, 0.80, f't\n{t}', fontsize=fstick+4, color='grey', ha='right')

    for axrow, df in zip(axs[1:], dfs):
        for ax, trait, traitvmax in zip(axrow, traits, traitvmaxs):
            df_piv = pd.pivot_table(df.loc[df.Time == t], values=trait, index=[pivindex], columns=['ES']).sort_index(axis=0, ascending=False)
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
