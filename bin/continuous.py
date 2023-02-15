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

traits = ['a2Seenmean',
            'ChooseGrainmean',
            'MimicGrainmean',
            'wmean']
traitlabels = ['Effort to get $\it{B}$',
                'Sensitivity for\nchoosing partner',
                'Sensitivity for\nmimicking partner',
                'Fitness']
folders = ['given0', 'none', 'p', 'r', 'pr', 'p8r']

movie = False

numa2 = 64
filename = 'output'

fslabel = 32 # Label font size
fstick = 24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

dfs = []
for folder in folders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.csv'))),
                    ignore_index=True)
    df.ChooseGrainmean = 1.0 - df.ChooseGrainmean
    df.MimicGrainmean = 1.0 - df.MimicGrainmean
    dfs.append(df)

df = dfs[1]
ts = df.Time.unique()
if movie:
    frames = []
else:
    ts = [ts[-1]]

given = pd.unique(df.Given)
if given > 0.9999999:
    given = 0.9999999
alphas = np.sort(pd.unique(df.alpha))[::-1]
rowindex = 'alpha'
logess = np.sort(pd.unique(df.logES))
rhos = 1.0 - 1.0/pow(2.0, logess)
nr = len(alphas)
nc = len(logess)
X, Y = np.meshgrid(np.linspace(0.0, mymodule.a2max, num=numa2),
                    np.linspace(mymodule.a2max, 0.0, num=numa2))
X = np.tile(A=X, reps=[nr, nc])
Y = np.tile(A=Y, reps=[nr, nc])
H = np.copy(Y)
Y[(X > Y)] = X[(X > Y)] 
X[(X > H)] = H[(X > H)] 
RRR, AAA = np.meshgrid(np.repeat(rhos, numa2),
                        np.repeat(alphas, numa2))
R = mymodule.fitness(Y, Y, given, AAA, RRR)
P = mymodule.fitness(X, X, given, AAA, RRR)
T = mymodule.fitness(Y, X, given, AAA, RRR)
S = mymodule.fitness(X, Y, given, AAA, RRR)

xmin = logess[0]
xmax = logess[-1]
xlabel = 'Substitutability of $\it{B}$'
ymin = alphas[-1]
ymax = alphas[0]
ylabel = 'Value of $\it{B}$'

traitvmaxs = [mymodule.a2max,
                mymodule.a2max,
                mymodule.a2max,
                mymodule.fitness(np.array([mymodule.a2max]),
                                    np.array([mymodule.a2max]),
                                    np.array([0.0]),
                                    np.array([0.9]),
                                    np.array([5.0]))]
xticklabels = [round(xmin),
                round((xmin + xmax)/2),
                round(xmax)]
yticklabels = [round(ymin, 1),
                round((ymin + ymax)/2, 1),
                round(ymax, 1)]
extent = 0, nc, 0, nr
extenta2 = 0, nc*numa2, 0, nr*numa2

fig, axs = plt.subplots(nrows=len(folders)+1,
                        ncols=len(traits),
                        figsize=(6*len(traits), 6*(len(folders)+1)))
#fig.delaxes(axs[0, 1])
#fig.delaxes(axs[0, 2])
#fig.delaxes(axs[0, 3])
fig.supxlabel(xlabel, x=0.513, y=0.06, fontsize=fslabel*1.5)
fig.supylabel(ylabel, x=0.05, y=0.493, fontsize=fslabel*1.5)

letter = ord('b')
for axrow in axs:
    for ax in axrow:
        if ax.get_subplotspec().is_first_row():
            ax.text(0,
                    nr*numa2*1.035,
                    'a',
                    fontsize=fslabel,
                    weight='bold')
            ax.set(xticks=[0, nc*numa2/2, nc*numa2],
                    yticks=[0, nr*numa2/2, nr*numa2],
                    xticklabels=[],
                    yticklabels=[])
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        else:
            ax.set(xticks=[0, nc/2, nc],
                            yticks=[0, nr/2, nr],
                            xticklabels=[],
                            yticklabels=[])
            if letter<= ord('z'):
                ax.text(0, 
                        nr*1.035,
                        chr(letter),
                        fontsize=fslabel,
                        weight='bold')
            else:
                ax.text(0,
                        nr*1.035,
                        'a' + chr(letter - 26),
                        fontsize=fslabel,
                        weight='bold')
            if ax.get_subplotspec().is_first_col():
                ax.set_yticklabels(yticklabels, fontsize=fstick) 
            letter += 1
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

Z = np.full([nr*numa2, nc*numa2, 4], mymodule.colormap['white'])
mymodule.snowdriftcolors(T, R, P, S, Z)
ax = axs[0, 0]
ax.imshow(Z, extent=extenta2)
ax.set_title('Games', pad=50.0, fontsize=fslabel)

Z = np.full([nr*numa2, nc*numa2], -1.0)
mask = mymodule.snowdrift(T, R, P, S) | mymodule.prisoner(T, R, P, S)
Z[mask] = (R[mask] > P[mask])
Z[mask] = (S[mask] - P[mask])
ax = axs[0, 1]
ax.imshow(Z, extent=extenta2, vmin=-1, vmax=0.05)
ax.set_title('Games', pad=50.0, fontsize=fslabel)

Z = np.full([nr*numa2, nc*numa2], -1.0)
mask = mymodule.prisoner(T, R, P, S) 
Z[mask] = (T[mask] + S[mask] - 2.0*R[mask])
ax = axs[0, 2]
ax.imshow(Z, extent=extenta2, vmin=-1, vmax=0.1)
ax.set_title('Games', pad=50.0, fontsize=fslabel)

for t in ts:
    for axrow, df in zip(axs[1:], dfs):
        for ax, trait, traitvmax in zip(axrow, traits, traitvmaxs):
            Z = pd.pivot_table(df.loc[df.Time == t],
                                values=trait,
                                index=[rowindex],
                                columns=['logES']).sort_index(axis=0,
                                                            ascending=False)
            ax.imshow(Z,
                    extent=extent,
                    cmap='viridis',
                    vmin=0,
                    vmax=traitvmax)
    if movie:
        text = fig.text(0.90,
                        0.93,
                        f't\n{t}',
                        fontsize=fstick+4,
                        color='grey',
                        ha='right')
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
