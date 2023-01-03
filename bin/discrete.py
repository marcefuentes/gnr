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

traits = ['a2Seenmean', 'ChooseGrainmean', 'MimicGrainmean', 'wmean']
traitlabels = ['Effort to get $\it{B}$', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner', 'Fitness']
folders = ['none', 'p', 'r', 'pr', 'p8r', 'given0']

movie = False

numa2 = 2
filename = 'output'

fslabel = 32 # Label font size
fstick = 24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

dfs = []
for folder in folders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.csv'))), ignore_index=True)
    df.ChooseGrainmean = 1.0 - df.ChooseGrainmean
    df.MimicGrainmean = 1.0 - df.MimicGrainmean
    df['help'] = df.a2Seenmean*mymodule.R2*df.Given
    dfs.append(df)

df = dfs[0]
ts = df.Time.unique()
if movie:
    frames = []
else:
    ts = [ts[-1]]
givens = np.sort(pd.unique(df.Given))[::-1]
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
    ymin = givens[-1]
    ymax = givens[0]
    ylabel = 'Partner\'s share of $\it{B}$'
    pivindex = 'Given'
else:
    nr = len(alphas)
    GG = np.full([nr, nc], givens[0])
if len(alphas) > 1:
    nr = len(alphas)
    RR, AA = np.meshgrid(rhos, alphas)
    ymin = alphas[-1]
    ymax = alphas[0]
    ylabel = 'Value of $\it{B}$'
    pivindex = 'alpha'
else:
    nr = len(givens)
    AA = np.full([nr, nc], alphas[0])

traitvmaxs = [mymodule.a2max, mymodule.a2max, mymodule.a2max, mymodule.fitness(np.array([mymodule.a2max]), np.array([mymodule.a2max]), np.array([0.0]), np.array([0.9]), np.array([5.0]))]
xticklabels = [round(xmin), round((xmin + xmax)/2), round(xmax)]
yticklabels = [round(ymin, 1), round((ymin + ymax)/2, 1), round(ymax, 1)]
extent = 0, nc, 0, nr
extenta2 = 0, nc, 0, nr*numa2

zeros = np.zeros([nr, nc])
a20 = np.copy(zeros)
a21 = a20 + mymodule.a2max/2.0
a22 = a20 + mymodule.a2max
a2eq = np.copy(zeros)
weq = np.copy(zeros)
xeq = np.copy(zeros)
Z0 = np.full([nr, nc, 4], mymodule.default)
Z1 = np.full([nr, nc, 4], mymodule.default)

w00 = mymodule.fitness(a20, a20, GG, AA, RR)
w01 = mymodule.fitness(a20, a21, GG, AA, RR)
w02 = mymodule.fitness(a20, a22, GG, AA, RR)
w10 = mymodule.fitness(a21, a20, GG, AA, RR)
w11 = mymodule.fitness(a21, a21, GG, AA, RR)
w12 = mymodule.fitness(a21, a22, GG, AA, RR)
w20 = mymodule.fitness(a22, a20, GG, AA, RR)
w21 = mymodule.fitness(a22, a21, GG, AA, RR)
w22 = mymodule.fitness(a22, a22, GG, AA, RR)

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

fig, axs = plt.subplots(nrows=len(folders)+1, ncols=len(traits), figsize=(6*len(traits), 6*(len(folders)+1)))
fig.delaxes(axs[0, 1])
fig.delaxes(axs[0, 2])
fig.delaxes(axs[0, 3])
fig.supxlabel(xlabel, x=0.513, y=0.06, fontsize=fslabel*1.5)
fig.supylabel(ylabel, x=0.05, y=0.493, fontsize=fslabel*1.5, ha='center')

letter = ord('b')
for axrow in axs:
    for ax in axrow:
        if ax.get_subplotspec().is_first_row():
            ax.set(xticks=[0, nc/2, nc], yticks=[0, nr*numa2/2, nr*numa2], xticklabels=[], yticklabels=[])
            ax.set_title('Game types', pad=50.0, fontsize=fslabel)
            ax.text(0, nr*numa2*1.035, 'a', fontsize=fslabel, weight='bold')
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

axs[0, 0].imshow(Z, extent=extenta2, aspect=1.0/numa2)

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
