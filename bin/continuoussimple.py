#! /usr/bin/env python

from glob import glob
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import pandas as pd
import time

start_time = time.perf_counter ()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

traits = ['ChooseGrainmean',
            'MimicGrainmean']
titles = ['Games',
                'Sensitivity for\nchoosing partner',
                'Sensitivity for\nmimicking partner']
folders = ['given95', 'given50']
givens = [0.95, 0.5]
subfolder = 'pr'

numa2 = 64

movie = False

fslarge = 32 # Label font size
fssmall = 24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

dfs = []
for folder in folders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, subfolder, '*.csv'))),
                    ignore_index=True)
    df.ChooseGrainmean = 1.0 - df.ChooseGrainmean
    df.MimicGrainmean = 1.0 - df.MimicGrainmean
    dfs.append(df)

df = dfs[0]
ts = df.Time.unique()
if movie:
    frames = []
else:
    ts = [ts[-1]]

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
R = mymodule.fitness(Y, Y, 0.0, AAA, RRR)
P = mymodule.fitness(X, X, 0.0, AAA, RRR)

xmin = logess[0]
xmax = logess[-1]
xlabel = 'Substitutability of $\it{B}$'
ymin = alphas[-1]
ymax = alphas[0]
ylabel = 'Value of $\it{B}$'

traitvmaxs = [mymodule.a2max,
                mymodule.a2max]
xticklabels = [round(xmin),
                round((xmin + xmax)/2),
                round(xmax)]
yticklabels = [round(ymin, 1),
                round((ymin + ymax)/2, 1),
                round(ymax, 1)]
extent = 0, nc, 0, nr
extenta2 = 0, nc*numa2, 0, nr*numa2

fig, axs = plt.subplots(nrows=len(givens),
                        ncols=len(titles),
                        figsize=(6*len(titles), 6*(len(givens))))
fig.supxlabel(xlabel, x=0.513, y=0.01, fontsize=fslarge*1.2)
fig.supylabel(ylabel, x=0.03, y=0.493, fontsize=fslarge*1.2)

letter = ord('a')
for axrow in axs:
    for ax, title in zip(axrow, titles):
        if title == 'Games':
            ax.text(0, 
                    nr*numa2*1.035,
                    chr(letter),
                    fontsize=fslarge*0.8,
                    weight='bold')
        else:
            ax.text(0, 
                    nr*1.035,
                    chr(letter),
                    fontsize=fslarge*0.8,
                    weight='bold')
        letter += 1
        if ax.get_subplotspec().is_first_col():
            ax.set(xticks=[0, nc*numa2/2, nc*numa2],
                    yticks=[0, nr*numa2/2, nr*numa2],
                    xticklabels=[])
            ax.set_yticklabels(yticklabels, fontsize=fssmall) 
        else:
            ax.set(xticks=[0, nc/2, nc],
                            yticks=[0, nr/2, nr],
                            xticklabels=[],
                            yticklabels=[])
        if ax.get_subplotspec().is_first_row():
            ax.set_title(title, pad=40.0, fontsize=fslarge*0.9)
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fssmall)

if givens[0] > 0.9999999:
    given = 0.9999999
for axrow, given in zip(axs, givens):

    T = mymodule.fitness(Y, X, given, AAA, RRR)
    S = mymodule.fitness(X, Y, given, AAA, RRR)
    Z = np.full([nr*numa2, nc*numa2, 4], mymodule.colormap['white'])
    mymodule.gamecolors(T, R, P, S, Z)
    axrow[0].imshow(Z, extent=extenta2)

for t in ts:
    for axrow, df in zip(axs, dfs):
        for ax, trait, traitvmax in zip(axrow[1:], traits, traitvmaxs):
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
                        fontsize=fssmall+4,
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
