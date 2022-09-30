#! /usr/bin/env python

from glob import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time

start_time = time.perf_counter ()

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

letters = [['a', 'b', 'c', 'd', 'e'],
            ['f', 'g', 'h', 'i', 'j'],
            ['k', 'l', 'm', 'n', 'o'],
            ['p', 'q', 'r', 's', 't'],
            ['u', 'v', 'w', 'x', 'y'],
            ['z', 'aa', 'ab', 'ac', 'ad']]

traits = ['a2Seenmedian', 'helpmedian', 'wmedian', 'ChooseGrainmedian', 'MimicGrainmedian']
traitlabels = ['Effort to get $\it{A}$', 'Help', 'Fitness', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner']
traitvmaxs = [0.5, 1.0, 1.0, 1.0, 1.0]
folders = ['none', 'p', 'r', 'pr', 'p8r']

alpha = 0.5
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

num = 2000
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

def fitness(A, Apartner):
    q1 = R1*(a2max - A)/b
    q2 = R2*(A*(1.0 - Ygivens) + Apartner*Ygivens)
    w = np.where(Xrhos == 0.0, pow(q1, alpha)*pow(q2, 1.0 - alpha), pow(alpha*pow(q1, Xrhos) + (1.0 - alpha)*pow(q2, Xrhos), 1.0/Xrhos)) 
    return w

def a2eq(given):
    if (R == 1.0) and (b == 1.0) and (alpha == 0.5):
        T = 1.0 - given
        Q = pow(T, 1.0/(Xrhos - 1.0))
    else:
        T = b*R*(1.0 - given)
        Q = R*pow(T*(1.0 - alpha)/alpha, 1.0/(Xrhos - 1.0))
    a2 = a2max/(1.0 + Q*b)
    return a2
    
R = R2/R1
b = a2max/a1max
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
givens = np.linspace(maxgiven, mingiven, num=num)
givens[0] = 0.99999
Xrhos, Ygivens = np.meshgrid(rhos, givens)
Aeq = a2eq(Ygivens)
Feq = fitness(Aeq, Aeq)
Zs = [Aeq, Aeq*R2*Ygivens, Feq, np.ones([num, num])*0.1, np.ones([num, num])*0.1]

fslabel=36 # Label font size
fstick=24 # Tick font size

fig, axs = plt.subplots(nrows=len(folders)+1, ncols=len(traits), figsize=(6*len(traits), 6*(len(folders)+1))) # constrained_layout=False, squeeze=False

fig.supxlabel('Substitutability of $\it{A}$', x=0.513, y=0.05, fontsize=fslabel*1.25)
fig.supylabel('Partner\'s share of $\it{A}$', x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

extent = 0, num, 0, num
for ax, Z, traitvmax, letter, traitlabel in zip(axs[0], Zs, traitvmaxs, letters[0], traitlabels):
    ax.imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)
    ax.set_xticks([0, num/2.0, num])
    ax.set_yticks([0, num/2.0, num])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if traitlabel == traitlabels[0]:
        ax.set_yticklabels([round(mingiven, 1), round((mingiven + maxgiven)/2, 1), round(maxgiven, 1)], fontsize=fstick) 
    ax.text(0, num*1.035, letter, fontsize=fslabel, weight='bold')

dfs = {}
for folder in folders:
    dfs[folder] = pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.csv'))), ignore_index=True)

t = dfs[folders[0]].Time.iat[-1]

dfts = {}
for folder in folders:
    dfts[folder] = dfs[folder].loc[dfs[folder]['Time'] == t].copy()
    dfts[folder].sort_values(by=['ES', 'Given'], inplace=True)
nr = len(pd.unique(dfts[folders[0]]['Given']))
nc = len(pd.unique(dfts[folders[0]]['ES']))

extent = 0, nr, 0, nc
for axrow, folder, letterrow in zip(axs[1:], folders, letters[1:]):
    for ax, trait, traitvmax, letter, traitlabel in zip(axrow, traits, traitvmaxs, letterrow, traitlabels):
        df = dfts[folder]
        if 'Sensitivity' in traitlabel: df[trait] = 1.0 - df[trait]
        df_piv = pd.pivot_table(df, values=trait, index=['Given'], columns=['ES']).sort_index(axis=0, ascending=False)
        ax.imshow(df_piv, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)
        ax.set_xticks([0, nc/2.0, nc])
        ax.set_yticks([0, nr/2.0, nr])
        ax.set_xticklabels([])
        if folder == folders[-1]:
            ax.set_xticklabels([round(minlog_es), round((minlog_es + maxlog_es)/2), round(maxlog_es)], fontsize=fstick)
        ax.set_yticklabels([])
        if traitlabel == traitlabels[0]:
            ax.set_yticklabels([round(mingiven, 1), round((mingiven + maxgiven)/2, 1), round(maxgiven, 1)], fontsize=fstick) 
        ax.text(0, nr*1.035, letter, fontsize=fslabel, weight='bold')

plt.savefig('coop.png', transparent=False)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
