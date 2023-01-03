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
datafolderp = 'discrete/p'
datafolderr = 'discrete/r'
movie = False

numa2 = 2
filename = 'discrete'

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

RR, GG = np.meshgrid(rhos, givens)

traitvmaxs = [mymodule.a2max, mymodule.a2max]
xticklabels = [round(xmin), round((xmin + xmax)/2), round(xmax)]
yticklabels = [round(ymin, 1), round((ymin + ymax)/2, 1), round(ymax, 1)]
extent = 0, nc, 0, nr
extenta2 = 0, nc, 0, nr*numa2

zeros = np.zeros([nr, nc])

fig, axs = plt.subplots(nrows=len(alphafolders), ncols=len(traitlabels), figsize=(6*len(traitlabels), 6*len(alphafolders)))
fig.supxlabel(xlabel, x=0.512, y=0.03, fontsize=fslabel*1.25)
fig.supylabel(ylabel, x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

letter = ord('a')
for axrow in axs:
    for ax, traitlabel in zip(axrow, traitlabels):
        if ax.get_subplotspec().is_first_row():
            ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)
        if ax.get_subplotspec().is_first_col():
            ax.set(xticks=[0, nc/2, nc], yticks=[0, nr*numa2/2, nr*numa2], xticklabels=[], yticklabels=[])
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
            ax.text(0, nr*numa2*1.035, chr(letter), fontsize=fslabel, weight='bold')
        else:
            ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
            ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
        letter += 1

for axrow, alpha in zip(axs, alphas):
    AA = np.full([nr, nc], alpha)
    a20 = np.copy(zeros)
    a21 = a20 + mymodule.a2max/2.0
    a22 = a20 + mymodule.a2max

    w00 = mymodule.fitness(a20, a20, GG, AA, RR)
    w01 = mymodule.fitness(a20, a21, GG, AA, RR)
    w02 = mymodule.fitness(a20, a22, GG, AA, RR)
    w10 = mymodule.fitness(a21, a20, GG, AA, RR)
    w11 = mymodule.fitness(a21, a21, GG, AA, RR)
    w12 = mymodule.fitness(a21, a22, GG, AA, RR)
    w20 = mymodule.fitness(a22, a20, GG, AA, RR)
    w21 = mymodule.fitness(a22, a21, GG, AA, RR)
    w22 = mymodule.fitness(a22, a22, GG, AA, RR)

    a2social = np.copy(zeros)
    mask = (w11 > w00)
    a2social[mask] = a21[mask]
    mask = (w22 > w11)
    a2social[mask] = a22[mask]
    wsocial = mymodule.fitness(a2social, a2social, GG, AA, RR)

    a2eq = np.copy(zeros)
    weq = np.copy(zeros)
    xeq = np.copy(zeros)
    Z0 = np.full([nr, nc, 4], mymodule.default)
    Z1 = np.full([nr, nc, 4], mymodule.default)

    mask0 = (w00 > w11)
    T = np.copy(w01)
    R = np.copy(w00)
    P = np.copy(w11)
    S = np.copy(w10)
    mymodule.gametypes(mask0, T, R, P, S, a20, a21, Z0, a2eq, xeq, weq)

    mask0 = (w00 < w11)
    T = np.copy(w10)
    R = np.copy(w11)
    P = np.copy(w00)
    S = np.copy(w01)
    mymodule.gametypes(mask0, T, R, P, S, a21, a20, Z0, a2eq, xeq, weq)

    mask0 = (w11 > w22)
    T = np.copy(w12)
    R = np.copy(w11)
    P = np.copy(w22)
    S = np.copy(w21)
    mymodule.gametypes(mask0, T, R, P, S, a21, a22, Z1, a2eq, xeq, weq)

    mask0 = (w11 < w22)
    T = np.copy(w21)
    R = np.copy(w22)
    P = np.copy(w11)
    S = np.copy(w12)
    mymodule.gametypes(mask0, T, R, P, S, a22, a21, Z1, a2eq, xeq, weq)

    Z = np.full([nr*numa2, nc, 4], mymodule.default)
    Z[::2,:] = Z1
    Z[1::2,:] = Z0

    axrow[0].imshow(Z, extent=extenta2, aspect=1.0/numa2)

for t in ts:

    if movie:
        text = fig.text(0.90, 0.90, f't\n{t}', fontsize=fstick+4, color='grey', ha='right')

    for axrow, dfp, dfr in zip(axs, dfps, dfrs):
        df = dfp.loc[dfp.Time == t].copy()
        df[traits[0]] = 1.0 - df[traits[0]]
        df_piv = pd.pivot_table(df, values=traits[0], index=[pivindex], columns=['logES']).sort_index(axis=0, ascending=False)
        axrow[1].imshow(df_piv, extent=extent, cmap='magma', vmin=0, vmax=traitvmaxs[0])
        df = dfr.loc[dfr.Time == t].copy()
        df[traits[1]] = 1.0 - df[traits[1]]
        df_piv = pd.pivot_table(df, values=traits[1], index=[pivindex], columns=['logES']).sort_index(axis=0, ascending=False)
        axrow[2].imshow(df_piv, extent=extent, cmap='magma', vmin=0, vmax=traitvmaxs[1])

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
