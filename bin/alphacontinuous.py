#! /usr/bin/env python

from glob import glob
import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import pandas as pd
import time

start_time = time.perf_counter ()

traits = ['ChooseGrainmean', 'MimicGrainmean']
traitlabels = ['Game types', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner']
alphafolders = ['alpha75', 'alpha50', 'alpha25']
alphas = [0.75, 0.50, 0.25]
datafolderp = 'continuous/p'
datafolderr = 'continuous/r'
movie = False

numa2 = 128
filename = 'continuous'

fslabel = 32 # Label font size
fstick = 24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

dfps = []
dfrs = []
for alphafolder in alphafolders:
    dfp = pd.concat(map(pd.read_csv, glob(os.path.join(alphafolder, datafolderp, '*.csv'))), ignore_index=True)
    dfr = pd.concat(map(pd.read_csv, glob(os.path.join(alphafolder, datafolderr, '*.csv'))), ignore_index=True)
    dfps.append(dfp)
    dfrs.append(dfr)

df = dfps[0]
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
xmin = logess[0]
xmax = logess[-1]
ymin = givens[-1]
ymax = givens[0]
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Partner\'s share of $\it{B}$'
pivindex = 'Given'

RRR, GGG = np.meshgrid(np.repeat(rhos, numa2), np.repeat(givens, numa2))

traitvmaxs = [mymodule.a2max, mymodule.a2max]
xticklabels = [round(xmin), round((xmin + xmax)/2), round(xmax)]
yticklabels = [round(ymin, 1), round((ymin + ymax)/2, 1), round(ymax, 1)]
extent = 0, nc, 0, nr
extenta2 = 0, nc*numa2, 0, nr*numa2

X, Y = np.meshgrid(np.linspace(0.0, mymodule.a2max, num=numa2), np.linspace(mymodule.a2max, 0.0, num=numa2))
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
            ax.set(xticks=[0, nc*numa2/2, nc*numa2], yticks=[0, nr*numa2/2, nr*numa2], xticklabels=[]) 
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
        else:
            ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
            ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
        letter += 1

for axrow, alpha in zip(axs, alphas):
    AAA = np.full([nr*numa2, nc*numa2], alpha)
    T = mymodule.fitness(Y, X, GGG, AAA, RRR)
    R = mymodule.fitness(Y, Y, GGG, AAA, RRR)
    P = mymodule.fitness(X, X, GGG, AAA, RRR)
    S = mymodule.fitness(X, Y, GGG, AAA, RRR)
    mask = (R < P)
    H = T[mask]
    T[mask] = S[mask]
    S[mask] = H
    H = R[mask]
    R[mask] = P[mask]
    P[mask] = H
    Z = np.full((nr*numa2, nc*numa2, 4), mymodule.colormap['default'])
    Z[(T < R) & (P < S)] = mymodule.colormap['nodilemma']
    Z[(T < R) & (P < S) & (2.0*R <= T + S)] = mymodule.colormap['nodilemmaRS']
    Z[(T > R) & (P > S)] = mymodule.colormap['prisoner']
    Z[(T > R) & (P > S) & (2.0*R <= T + S)] = mymodule.colormap['prisonerRS']
    Z[(T >= R) & (P <= S)] = mymodule.colormap['snowdrift']
    Z[(T >= R) & (P <= S) & (2.0*R <= T + S)] = mymodule.colormap['snowdriftRS']
    Z[R == P] = mymodule.colormap['nodilemma']

    axrow[0].imshow(Z, extent=extenta2)

for t in ts:

    if movie:
        text = fig.text(0.90, 0.90, f't\n{t}', fontsize=fstick+4, color='grey', ha='right')

    for axrow, dfp, dfr in zip(axs, dfps, dfrs):
        df = dfp.loc[dfp.Time == t].copy()
        df[traits[0]] = 1.0 - df[traits[0]]
        Z = pd.pivot_table(df, values=traits[0], index=[pivindex], columns=['logES']).sort_index(axis=0, ascending=False)
        axrow[1].imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=traitvmaxs[0])
        df = dfr.loc[dfr.Time == t].copy()
        df[traits[1]] = 1.0 - df[traits[1]]
        Z = pd.pivot_table(df, values=traits[1], index=[pivindex], columns=['logES']).sort_index(axis=0, ascending=False)
        axrow[2].imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=traitvmaxs[1])

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
