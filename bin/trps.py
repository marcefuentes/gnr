#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = f'trps{sys.argv[1]}.png'

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
log_ess = np.linspace(-5, 5, num=11) if sys.argv[1] == 'ces' else np.linspace(-10, 0, num=11)
ess = pow(2, log_ess)
givens = np.linspace(1.0, 0.0, num=11)
titles = []
for aC, aD in zip(aCs, aDs): 
    titles.append('$\it{a}$ = {' + str(aD) + ', ' + str(aC) + '}')
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

fig, axs = plt.subplots(nrows=1, ncols=len(aCs), sharey=True, figsize=(5.0*len(aCs) + 1.0, height), constrained_layout=False, squeeze=False)
axs = axs.flatten()
fig.supxlabel('Substitutability of $\it{A}$', fontsize=fslabel)
fig.supylabel('Partner\'s share of $\it{A}$', fontsize=fslabel)

for ax, title, letter, aC, aD in zip(axs, titles, letters, aCs, aDs):
    q1C = R1*(1.0-aC);
    q1D = R1*(1.0-aD);
    for es, log_es in zip(ess, log_ess):
        rho = 1.0 - 1.0/es if sys.argv[1] == 'ces' else es
        R = wC1 = fitness(q1C, R2*aC, rho)
        P = wD0 = fitness(q1D, R2*aD, rho)
        for given in givens:
            x = []
            y = []
            T = wD1 = fitness(q1D, R2*(aD*(1.0-given) + aC*given), rho)
            S = wC0 = fitness(q1C, R2*(aC*(1.0-given) + aD*given), rho)
            trpss = np.linspace(log_es-0.3, log_es+0.3, num=4)
            [x.append(pow(2, trps)) for trps in trpss]
            y.append(given + (T-0.8)/6.0)
            y.append(given + (R-0.8)/6.0)
            y.append(given + (P-0.8)/6.0)
            y.append(given + (S-0.8)/6.0)
            shift = shiftconstant*(S-P)
            if (T<R) and (P<S):
                # rgb = (0.5-shift, 0.5+shift, 0.5+shift)
                rgb = (0.5, 0.0, 1.0)
            else:
                if (T>R) and (P<S):
                    rgb = (0.0, 1.0, 0.0) # Snowdrift
                else:
                    rgb = (1.0, 0.0, 0.5)
                    # rgb = (0.5-shift, 0.5+shift, 0.5+shift)
            ax.plot(x, y, color=rgb, marker='o', markerfacecolor='white', linewidth=1.0, markersize=3)
            ax.set_title(title, pad=10.0, fontsize=fstitle)
            ax.tick_params(axis='x', labelsize=fstick)
            ax.tick_params(axis='y', labelsize=fstick)
            ax.set_xlim(pow(2, -5.5), pow(2, 5.5))
            ax.set_ylim(-0.08, 1.08)
            ax.set_xscale('log', base=2)
            ax.set_box_aspect(1)
            ax.text(0.01, 1.1, letter, fontsize=18, weight='bold')

fig.savefig(filename, transparent=False)
plt.close()

