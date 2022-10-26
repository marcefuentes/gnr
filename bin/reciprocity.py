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
traitvmaxs = [0.5, 1.0, 1.0, 1.0, 0.25]
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
r = 1.0
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
T = fitness(a2Cs, a2Ds)
R = fitness(a2Cs, a2Cs)
P = fitness(a2Ds, a2Ds)
S = fitness(a2Ds, a2Cs)
x = np.full([nc, nr], 0.0)
t = np.full([nc, nr], 0.0)
#mask = (T < R) & (P < S)
#x[mask] = 1.0
#t[mask] = 0.0
denominator = P - T - P*r + T*r
mask = denominator != 0.0
x[mask] = (P[mask] - T[mask] - P[mask]*r + R[mask]*r)/denominator[mask]
t[mask] = (T[mask] - R[mask]*r)/denominator[mask]
mask = x < 0.0
x[mask] = 0.0
mask = t < 0.0
t[mask] = 0.0
mask = x > 1.0
x[mask] = 1.0
mask = t > 1.0
t[mask] = 1.0
y = 1.0 - x - t
a2 = (x + t*(x + t) + t*y/r)*a2max/2.0
helps = a2*R2*GG 
wA = R*(x + t) + S*y
wT = R*(x + t) + S*y/r + P*y*(r - 1.0)/r
wB = T*x + T*t/r + P*t*(r - 1.0)/r + P*y
w = wA*x + wT*t + wB*y 
Zs = [a2, helps, w, np.ones([nc, nr])*0.06, t]

fslabel=36 # Label font size
fstick=24 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig, axs = plt.subplots(nrows=1, ncols=len(traits), figsize=(6*len(traits), 6))

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
plt.savefig('reciprocity.png', transparent=False)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
