#! /usr/bin/env python

from glob import glob
from math import log
import numpy as np
import os
import time
import matplotlib.pyplot as plt
import pandas as pd

start_time = time.perf_counter ()

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

letters = [['a', 'b', 'c', 'd', 'e'],
            ['f', 'g', 'h', 'i', 'j'],
            ['k', 'l', 'm', 'n', 'o'],
            ['p', 'q', 'r', 's', 't'],
            ['u', 'v', 'w', 'x', 'y']]

folders = ['none', 'p', 'r', 'pr']
traits = ['a2Seenmedian', 'helpmedian', 'wmedian', 'ChooseGrainmedian', 'MimicGrainmedian']
traitlabels = ['Effort to get $\it{A}$', 'Help', 'Fitness', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner']
traitvmaxs = [0.5, 1.0, 1.0, 1.0, 1.0]

alpha = 0.5
R1 = 2.0
R2 = 2.0
R = R2/R1
a1max = 1.0
a2max = 1.0
b = a2max/a1max

num = 2000
log_ess = np.linspace(-5, 5, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
givens = np.linspace(1.0, 0.0, num=num)
givens[0] = 0.99999
Xrhos, Ygivens = np.meshgrid(rhos, givens)

def a2eq(X, Y):
    T = b*R*(1.0 - Y)
    Q = R*pow(T*(1.0 - alpha)/alpha, 1.0/(X - 1.0))
    a2 = 1.0/(1.0 + Q)
    return a2

def fitness(Z, X, Y):
    q1 = (1.0 - Z)*R1
    q2 = Z*R2*(1.0 - Y) + Z*R2*Y
    w = np.where(X == 0.0, pow(q1, alpha)*pow(q2, 1.0 - alpha), pow(alpha*pow(q1, X) + (1.0 - alpha)*pow(q2, X), 1.0/X)) 
    return w

Z0 = a2eq(Xrhos, Ygivens)
Zs = [Z0, Z0*R2*Ygivens, fitness(Z0, Xrhos, Ygivens), np.ones([num, num])*0.1, np.ones([num, num])*0.1]
Zmaxs = [0.5, 1.0, 1.0, 1.0, 1.0]

dfs = {}
for folder in folders:
    dfs[folder] = pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.csv'))), ignore_index=True)

t = dfs[folders[0]].Time.iat[-1]

dfts = {}
for folder in folders:
    dfts[folder] = dfs[folder].loc[dfs[folder]['Time'] == t].copy()
    dfts[folder].sort_values(by=['ES', 'Given'], inplace=True)

x = dfts[folders[0]]['ES']
x = x.apply(lambda i: log(i, 2))
y = dfts[folders[0]]['Given']

fslabel = 26 # Label font size
fstick = 18 # Tick font size

fig, axs = plt.subplots(nrows=5, ncols=5, figsize=(30, 30)) # constrained_layout=False, squeeze=False

fig.supxlabel('Substitutability of $\it{A}$', x=0.513, y=0.05, fontsize=fslabel*1.5)
fig.supylabel('Partner\'s share of $\it{A}$', x=0.05, y=0.493, fontsize=fslabel*1.5, ha='center')

extent = 0, num, 0, num
for ax, Z, Zmax, letter, traitlabel in zip(axs[0], Zs, Zmaxs, letters[0], traitlabels):
    ax.imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=Zmax)
    ax.set_title(traitlabel, pad=40.0, fontsize=fslabel)
    ax.set_xticks([0, num/2.0, num])
    ax.set_yticks([0, num/2.0, num])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if 'Effort' in traitlabel:
        ax.set_yticklabels([0.0, 0.5, 1.0], fontsize=fstick) 
    ax.text(0, 2070, letter, fontsize=fslabel, weight='bold')

for axrow, folder, letterrow in zip(axs[1:], folders, letters[1:]):
    for ax, trait, traitvmax, letter, traitlabel in zip(axrow, traits, traitvmaxs, letterrow, traitlabels):
        df = dfts[folder]
        s = df[trait]
        if 'Sensitivity' in traitlabel: s = 1 - s
        ax.scatter(x, y, c=s, cmap='magma', marker='s', vmin=0, vmax=traitvmax, s=700)
        ax.set_xlim(-5.5, 5.5)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xticks([-5, 0, 5])
        ax.set_yticks([0.0, 0.5, 1.0])
        ax.set_xticklabels([])
        if 'pr' in folder:
            ax.set_xticklabels([-5, 0, 5], fontsize=fstick)
        ax.set_yticklabels([])
        if 'Effort' in traitlabel:
            ax.set_yticklabels([0.0, 0.5, 1.0], fontsize=fstick) 
        ax.text(-5.5, 1.09, letter, fontsize=fslabel, weight='bold')
        ax.set_box_aspect(1)

plt.savefig('coop.png', transparent=False)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
