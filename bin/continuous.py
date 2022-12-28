#! /usr/bin/env python

from glob import glob
import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time

start_time = time.perf_counter ()

traits = ['a2Seenmean', 'ChooseGrainmean', 'MimicGrainmean', 'wmean']
traitlabels = ['Effort to get $\it{B}$', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner', 'Fitness']
folders = ['none', 'pr']

movie = False

numa2 = 64
filename = 'output'
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
logess = np.sort(pd.unique(df.logES))
rhos = 1.0 - 1.0/pow(2.0, logess)
alphas = np.sort(pd.unique(df.alpha))[::-1]
nc = len(rhos)
xmin = logess[0]
xmax = logess[-1]
xlabel = 'Substitutability of $\it{B}$'

if len(givens) > 1:
    nr = len(givens)
    RR, GG = np.meshgrid(rhos, givens)
    RRR, GGG = np.meshgrid(np.repeat(rhos, numa2), np.repeat(givens, numa2))
    ymin = givens[-1]
    ymax = givens[0]
    ylabel = 'Partner\'s share of $\it{B}$'
    pivindex = 'Given'
else:
    nr = len(alphas)
    GG = np.full([nc, nr], givens[0])
    GGG = np.full([nc*numa2, nr*numa2], givens[0])
if len(alphas) > 1:
    nr = len(alphas)
    RR, AA = np.meshgrid(rhos, alphas)
    RRR, AAA = np.meshgrid(np.repeat(rhos, numa2), np.repeat(alphas, numa2))
    ymin = alphas[-1]
    ymax = alphas[0]
    ylabel = 'Value of $\it{B}$'
    pivindex = 'alpha'
else:
    nr = len(givens)
    AA = np.full([nc, nr], alphas[0])
    AAA = np.full([nc*numa2, nr*numa2], alphas[0])

b = a2max/a1max
Rq = R2/R1
traitvmaxs = [a2max, a2max, a2max, fitness(np.array([a2max]), np.array([a2max]), np.array([0.0]), np.array([0.9]), np.array([5.0]))]
xticklabels = [round(xmin), round((xmin + xmax)/2), round(xmax)]
yticklabels = [round(ymin, 1), round((ymin + ymax)/2, 1), round(ymax, 1)]
extent = 0, nr, 0, nc
extenta2 = 0, nr*numa2, 0, nc*numa2
prisoner = [0.5, 0.0, 0.0, 1.0]
snowdrift = [0.0, 1.0, 1.0, 1.0]
nodilemma = [1.0, 1.0, 1.0, 1.0]
green = [0.0, 1.0, 0.0, 1.0]

MRT0 = b*Rq
Z = np.full((numa2*nr, numa2*nc, 4), green)
X, Y = np.meshgrid(np.linspace(0.0, a2max, num=numa2), np.linspace(a2max, 0.0, num=numa2))
X = np.tile(A=X, reps=[nr, nc])
Y = np.tile(A=Y, reps=[nr, nc])
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
Z[(T > R) & (P > S)] = prisoner
Z[(T >= R) & (P <= S) & (R != P)] = snowdrift
Z[((T < R) & (P < S)) | (R == P)] = nodilemma
#Z = np.tril(Z, k=-1)
#Z = np.ma.masked_where(Z == 0.0, Z)

MRT = MRT0*(1.0 - GG)
Q0 = Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2social = a2max/(1.0 + Q0*b)
Q = Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2eq = a2max/(1.0 + Q*b)
wsocial = fitness(a2social, a2social, GG, AA, RR)
weq = fitness(a2eq, a2eq, GG, AA, RR)
Mss = [[a2social, wsocial], [a2eq, weq]]

fig, axs = plt.subplots(nrows=len(folders)+1, ncols=len(traits), figsize=(6*len(traits), 6*(len(folders)+1)))
fig.delaxes(axs[0, 1])
fig.delaxes(axs[0, 2])
fig.delaxes(axs[0, 3])
fig.supxlabel(xlabel, x=0.513, y=0.03, fontsize=fslabel*1.25)
fig.supylabel(ylabel, x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

letter = ord('b')
for axrow in axs:
    for ax in axrow:
        if ax.get_subplotspec().is_first_row():
            ax.set(xticks=[0, numa2*nc/2, numa2*nc], yticks=[0, numa2*nr/2, numa2*nr], xticklabels=[], yticklabels=[])
            ax.set_title('Game types', pad=50.0, fontsize=fslabel)
            ax.text(0, numa2*nr*1.035, 'a', fontsize=fslabel, weight='bold')
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        else:
            ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
            ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
            letter += 1
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

axs[0, 0].imshow(Z, extent=extenta2)

for t in ts:

    if movie:
        text = fig.text(0.90, 0.90, f't\n{t}', fontsize=fstick+4, color='grey', ha='right')

    for axrow, df in zip(axs[1:], dfs):
        for ax, trait, traitvmax in zip(axrow, traits, traitvmaxs):
            df_piv = pd.pivot_table(df.loc[df.Time == t], values=trait, index=[pivindex], columns=['logES']).sort_index(axis=0, ascending=False)
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
