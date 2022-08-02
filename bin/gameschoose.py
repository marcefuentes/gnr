#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = f'wcwdchoose{sys.argv[1]}.png'

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

height = 6.0
fslabel = 26
fstitle= 24
fstick = 16    
shiftconstant = 5.0 # Micromutations
shiftconstant = 0.7 # Macromutations

# Parameters

alpha = 0.5
R1 = 2.0
R2 = 2.0

aCs = [0.3, 0.5]
aDs = [0.2, 0.4]
fx = 0.1
ess = np.linspace(-5, 5, num=11) if sys.argv[1] == 'ces' else np.linspace(-10, 0, num=11)
ess = pow(2, ess)
givens = np.linspace(1.0, 0.0, num=11)
titles = []
for aC, aD in zip(aCs, aDs): 
    titles.append('$\it{R}$ - $\it{P}$\nfor $\it{a}$ = {' + str(aD) + ', ' + str(aC) + '}')
letters = ['a', 'b', 'c']

def fitness(q1, q2, rho):
    if sys.argv[1] == 'q':
        w = 4.0*pow(q1, rho)/9.0 + 4.0*q2/9.0
    else:
        if rho == 0.0:
            w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
        else:
            w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

x = []
y = []
for es in ess:
    for given in givens:
        x.append(es)
        y.append(given)

sizes = [[], []]
colors = [[], []]
edgecolors = [[], []]
for aC, aD, size, color, edgecolor in zip(aCs, aDs, sizes, colors, edgecolors):
    q1C = R1*(1.0-aC);
    q1D = R1*(1.0-aD);
    for es in ess:
        rho = 1.0 - 1.0/es if sys.argv[1] == 'ces' else es
        R = wC1 = fitness(q1C, R2*aC, rho)
        P = wD0 = fitness(q1D, R2*aD, rho)
        for given in givens:
            T = wD1 = fitness(q1D, R2*(aD*(1.0-given) + aC*given), rho)
            S = wC0 = fitness(q1C, R2*(aC*(1.0-given) + aD*given), rho)
            wcwd0 = R*fx + S*(1.0-fx) - P*(1.0-fx) - T*fx
            wcwd = R - P
            rgb = (0.2, 1.0, 0.0) if (wcwd0 < 0.02) and (wcwd > 0) else (0.8, 1.0-wcwd0, 0.8)
            color.append(rgb)
            edgecolor.append(rgb)
            size.append(2000.0*abs(wcwd))

fig, axs = plt.subplots(nrows=1, ncols=len(aCs), sharey=True, figsize=(5.0*len(aCs) + 1.0, height), constrained_layout=False, squeeze=False)
axs = axs.flatten()
fig.supxlabel('Substitutability of $\it{A}$', fontsize=fslabel)
fig.supylabel('Partner\'s share of $\it{A}$', fontsize=fslabel)

for ax, size, color, edgecolor, title, letter in zip(axs, sizes, colors, edgecolors, titles, letters): 
    ax.scatter(x=x, y=y, s=size, color=color, ec=edgecolor)
    ax.set_title(title, pad=10.0, fontsize=fstitle)
    ax.tick_params(axis='x', labelsize=fstick)
    ax.tick_params(axis='y', labelsize=fstick)
    ax.set_xscale('log', base=2)
    ax.set_box_aspect(1)
    ax.text(0.01, 1.1, letter, fontsize=18, weight='bold')

fig.savefig(filename, transparent=False)
plt.close()

