#! /usr/bin/env python

from glob import glob
from math import log
import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time

start_time = time.perf_counter ()

movie = False

letters = [['a', 'b', 'c', 'd', 'e'],
            ['f', 'g', 'h', 'i', 'j'],
            ['k', 'l', 'm', 'n', 'o'],
            ['p', 'q', 'r', 's', 't'],
            ['u', 'v', 'w', 'x', 'y'],
            ['z', 'aa', 'ab', 'ac', 'ad']]

traits = ['a2Seenmean', 'help', 'wmean', 'ChooseGrainmean', 'MimicGrainmean']
traitlabels = ['Effort to get $\it{A}$', 'Help', 'Fitness', 'Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner']
traitvmaxs = [0.5, 1.0, 1.0, 1.0, 1.0]
folders = ['none', 'p', 'r', 'pr', 'p8r']

alpha = 0.5
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

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

r = 1.0/(1.0 - pow(1.0 - pow(2, -7), 2))
R = R2/R1
b = a2max/a1max
givens = np.linspace(1.0, 0.0, num=21)
givens[0] = 0.9999999
log_ess = np.linspace(-5.0, 5.0, num=21)
rhos = 1.0 - 1.0/pow(2, log_ess)
nr = len(givens)
nc = len(rhos)
RR, GG = np.meshgrid(rhos, givens)
a2Ds = np.full([nc, nr], 0.0)
a2Cs = np.full([nc, nr], a2max/2.0)
Ts = fitness(a2Cs, a2Ds)
Rs = fitness(a2Cs, a2Cs)
Ps = fitness(a2Ds, a2Ds)
Ss = fitness(a2Ds, a2Cs)
xeqss = np.full([nc, nr], 0.0)
mask = (Ts < Rs) & (Ps < Ss)
xeqss[mask] = 1.0
mask = (Ts > Rs) & (Ps < Ss) & (Rs - Ss - Ts + Ps != 0.0)
xeqss[mask] = (Ps[mask] - Ss[mask])/(Rs[mask] - Ss[mask] - Ts[mask] + Ps[mask])
weqss = (Ts + Ss)*xeqss*(1.0 - xeqss) + Rs*xeqss*xeqss + Ps*(1.0 - xeqss)*(1.0 - xeqss)
a2eqss = xeqss*a2max/2.0
helpeqss = a2eqss*R2*GG 
frTs = (2*Ps - Ss/r + Ps/r)/(Rs + Ps - Ss/r + Ts/r + Ps/r)
mask = frTs < 0.0
frTs[mask] = 0.0
Zs = [a2eqss, helpeqss, weqss, np.ones([nc, nr])*0.06, frTs]

fslabel=36 # Label font size
fstick=24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig, axs = plt.subplots(nrows=1, ncols=len(traits), figsize=(6*len(traits), 6*(len(folders)+1)))

fig.supxlabel('Substitutability of $\it{A}$', x=0.513, y=0.05, fontsize=fslabel*1.25)
fig.supylabel('Partner\'s share of $\it{A}$', x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

# All plots

extent = 0, nr, 0, nc
givens[0] = 1.0
minx = round(log_ess[0], 2)
maxx = round(log_ess[-1], 2)
miny = givens[-1]
maxy = givens[0]
xticklabels = [minx, round((minx + maxx)/2), maxx]
yticklabels = [miny, (miny + maxy)/2, maxy]

for ax in axs:
    ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
axs[0].set_yticklabels(yticklabels, fontsize=fstick) 

# Top row of plots

for ax, Z, traitvmax in zip(axs, Zs, traitvmaxs):
    ax.imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)
plt.savefig('alphar.png', transparent=False)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
