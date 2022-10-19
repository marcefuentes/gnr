#! /usr/bin/env python

from glob import glob
from math import log
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time

start_time = time.perf_counter ()

letters = [['a', 'b', 'c', 'd', 'e'],
            ['f', 'g', 'h', 'i', 'j'],
            ['k', 'l', 'm', 'n', 'o'],
            ['p', 'q', 'r', 's', 't'],
            ['u', 'v', 'w', 'x', 'y']]

traits = ['a2Seenmean', 'help', 'wmean']
traitlabels = ['Effort to get $\it{A}$', 'Help', 'Fitness']
traitvmaxs = [0.5, 1.0, 1.0]

a1max = 1.0
a2max = 1.0
num = 21         # Number of pixels
numalpha = 3    # Number of plot rows and columns
every = int(num/2)
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

def fitness(x, y):
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - GG) + x*R2*GG
    w = q1*q2
    mask = (w > 0.0) & (RR == 0.0)
    w[mask] = pow(q1[mask], alpha)*pow(q2[mask], 1.0 - alpha)
    mask = (w > 0.0) & (RR < 0.0)
    w[mask] = alpha*pow(q1[mask], RR[mask]) + (1.0 - alpha)*pow(q2[mask], RR[mask])
    mask = (w > 0.0) & (RR < 0.0)
    w[mask] = pow(w[mask], 1.0/RR[mask])
    mask = (RR > 0.0)
    w[mask] = pow(alpha*pow(q1[mask], RR[mask]) + (1.0 - alpha)*pow(q2[mask], RR[mask]), 1.0/RR[mask])
    return w

b = a2max/a1max
givens = np.linspace(maxgiven, mingiven, num=num)
givens[0] = 0.9999999
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
nr = len(givens)
nc = len(rhos)
RR, GG = np.meshgrid(rhos, givens)
alphas = np.linspace(1.0/(numalpha + 1), numalpha/(numalpha + 1), num=numalpha)
R2s = np.linspace(4.0/(numalpha + 1), numalpha*4.0/(numalpha + 1), num=numalpha)

fslabel=36 # Label font size
fstick=24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig, axs = plt.subplots(nrows=numalpha, ncols=numalpha, figsize=(6*numalpha, 6*numalpha))

fig.supxlabel('Substitutability of $\it{A}$', x=0.513, y=0.05, fontsize=fslabel*1.25)
fig.supylabel('Partner\'s share of $\it{A}$', x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

extent = 0, nr, 0, nc

for axrow, R2 in zip(axs, R2s):
    for ax, alpha in zip(axrow, alphas):
        R1 = 4.0 - R2
        R = R2/R1
        TT = b*R*(1.0 - GG)
        QQ = R*pow(TT*(1.0 - alpha)/alpha, 1.0/(RR - 1.0))
        a2eqss = a2max/(1.0 + QQ*b)
        #helpeqss = a2eqss*R2*GG
        #weqss = fitness(a2eqss, a2eqss)
        ax.imshow(a2eqss, extent=extent, cmap='magma', vmin=0, vmax=1.0)

# All plots

for axrow in axs:
    for ax in axrow:
        ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])

givens[0] = 1.0
minx = round(log_ess[0], 2)
maxx = round(log_ess[-1], 2)
miny = givens[-1]
maxy = givens[0]

xticklabels = [minx, round((minx + maxx)/2), maxx]
yticklabels = [miny, (miny + maxy)/2, maxy]
for axrow in axs:
    axrow[0].set_yticklabels(yticklabels, fontsize=fstick) 
for ax in axs[-1]:
    ax.set_xticklabels(xticklabels, fontsize=fstick)
for ax, alpha in zip(axs[-1], alphas):
    ax.set_xlabel(round(alpha, 2), fontsize=fslabel)
for axrow, R2 in zip(axs, R2s):
    axrow[0].set_ylabel(round(R2, 2), fontsize=fslabel)
for axrow, letterrow in zip(axs, letters):
    for ax, letter in zip(axrow, letterrow):
        ax.text(0, nr*1.035, letter, fontsize=fslabel, weight='bold')

plt.savefig('alpha.png', transparent=False)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
