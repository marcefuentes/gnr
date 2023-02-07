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

numa2 = 128
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
    #df['help'] = df.a2Seenmean*mymodule.R2*df.Given
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
    rowindex = 'Given'
else:
    nr = len(alphas)
    GG = np.full([nr, nc], givens[0])
    GGG = np.full([nr*numa2, nc*numa2], givens[0])
if len(alphas) > 1:
    nr = len(alphas)
    RR, AA = np.meshgrid(rhos, alphas)
    RRR, AAA = np.meshgrid(np.repeat(rhos, numa2), np.repeat(alphas, numa2))
    ymin = alphas[-1]
    ymax = alphas[0]
    ylabel = 'Value of $\it{B}$'
    rowindex = 'alpha'
else:
    nr = len(givens)
    AA = np.full([nr, nc], alphas[0])
    AAA = np.full([nr*numa2, nc*numa2], alphas[0])

traitvmaxs = [mymodule.a2max, mymodule.a2max, mymodule.a2max, mymodule.fitness(np.array([mymodule.a2max]), np.array([mymodule.a2max]), np.array([0.0]), np.array([0.9]), np.array([5.0]))]
xticklabels = [round(xmin), round((xmin + xmax)/2), round(xmax)]
yticklabels = [round(ymin, 1), round((ymin + ymax)/2, 1), round(ymax, 1)]
extent = 0, nc, 0, nr
extenta2 = 0, nc*numa2, 0, nr*numa2

MRT0 = mymodule.b*mymodule.Rq
X, Y = np.meshgrid(np.linspace(0.0, mymodule.a2max, num=numa2), np.linspace(mymodule.a2max, 0.0, num=numa2))
X = np.tile(A=X, reps=[nr, nc])
Y = np.tile(A=Y, reps=[nr, nc])
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

MRT = MRT0*(1.0 - GG)
Q0 = mymodule.Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
weq = mymodule.fitness(a2eq, a2eq, GG, AA, RR)

fig, axs = plt.subplots(nrows=len(folders)+1, ncols=len(traits), figsize=(6*len(traits), 6*(len(folders)+1)))
fig.delaxes(axs[0, 0])
fig.delaxes(axs[0, 3])
fig.supxlabel(xlabel, x=0.513, y=0.06, fontsize=fslabel*1.50)
fig.supylabel(ylabel, x=0.05, y=0.493, fontsize=fslabel*1.50, ha='center')

letter = ord('c')
for axrow in axs:
    for ax in axrow:
        if ax.get_subplotspec().is_first_row():
            ax.set(xticks=[0, nc*numa2/2, nc*numa2], yticks=[0, nr*numa2/2, nr*numa2], xticklabels=[], yticklabels=[])
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        else:
            ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
            if letter<= ord('z'):
                ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
            else:
                ax.text(0, nr*1.035, 'a' + chr(letter - 26), fontsize=fslabel, weight='bold')
            letter += 1
            if ax.get_subplotspec().is_first_col():
                ax.set_yticklabels(yticklabels, fontsize=fstick) 
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

#Z = np.full((nr*numa2, nc*numa2, 4), mymodule.colormap['nodilemma'])
#Z[(T < R) & (P < S) & (2.0*R <= T + S)] = mymodule.colormap['nodilemmaRS']
#axs[0, 0].imshow(Z, extent=extenta2)

Z = np.full((nr*numa2, nc*numa2, 4), mymodule.colormap['nodilemma'])
Z[(T >= R) & (P <= S)] = mymodule.colormap['snowdrift']
#Z[(T >= R) & (P <= S) & (2.0*R <= T + S)] = mymodule.colormap['snowdriftRS']
axs[0, 1].imshow(Z, extent=extenta2)
axs[0, 1].set_title('Snowdrift', pad=50.0, fontsize=fslabel)
axs[0, 1].text(0, nr*numa2*1.035, 'a', fontsize=fslabel, weight='bold')
axs[0, 1].set_yticklabels(yticklabels, fontsize=fstick) 

TS = np.zeros([nr*numa2, nc*numa2])
mask = (T > R) & (P > S)
TS[mask] =  (1.0 + T[mask] + S[mask] - 2.0*R[mask])/2.0 
axs[0, 2].imshow(TS, extent=extenta2, cmap='magma', vmin=0.0, vmax=0.7)
axs[0, 2].set_title('2$\it{R}$ - $\it{T}$ - $\it{S}$ under\nprisoner\'s dilemma', pad=50.0, fontsize=fslabel)
axs[0, 2].text(0, nr*numa2*1.035, 'b', fontsize=fslabel, weight='bold')

for t in ts:
    for axrow, df in zip(axs[1:], dfs):
        for ax, trait, traitvmax in zip(axrow, traits, traitvmaxs):
            Z = pd.pivot_table(df.loc[df.Time == t], values=trait, index=[rowindex], columns=['logES']).sort_index(axis=0, ascending=False)
            ax.imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)
    if movie:
        text = fig.text(0.90, 0.90, f't\n{t}', fontsize=fstick+4, color='grey', ha='right')
        plt.savefig('temp.png', transparent=False)
        text.remove()
        frames.append(iio.imread('temp.png'))
        os.remove('temp.png')
    else:
        plt.savefig(filename + '.png', transparent=False)

plt.close()

if movie:
    iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
