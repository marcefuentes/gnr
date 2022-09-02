#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = f'tr{sys.argv[1]}.png'

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

n = 16
aCs = np.linspace(0.5, 0.2, num=n)
#log_ess = np.linspace(-5, 5, num=11) if sys.argv[1] == 'ces' else np.linspace(-10, 0, num=11)
log_ess = np.linspace(-3, -3, num=1) if sys.argv[1] == 'ces' else np.linspace(-10, 0, num=11)
ess = pow(2, log_ess)
#givens = np.linspace(1.0, 0.0, num=11)
givens = [ 0.7 ]
titles = ['T - R', 'P - S']
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

fig, axs = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=(5.0*1 + 1.0, height), constrained_layout=False, squeeze=False)
axs = axs.flatten()
fig.supxlabel('Substitutability of $\it{A}$', fontsize=fslabel)
fig.supylabel('Partner\'s share of $\it{A}$', fontsize=fslabel)

for ax, title, letter in zip(axs, titles, letters):
    for es, log_es in zip(ess, log_ess):
        rho = 1.0 - 1.0/es if sys.argv[1] == 'ces' else es
        for given in givens:
            for aC in aCs:
                q1C = R1*(1.0-aC);
                R = wC1 = fitness(q1C, R2*aC, rho)
                aDs = np.linspace(aC-0.1, 0.1, num=int(aC*10) - 1)
                for aD in aDs:
                    q1D = R1*(1.0-aD);
                    P = wD0 = fitness(q1D, R2*aD, rho)
                    x = []
                    y = []
                    T = wD1 = fitness(q1D, R2*(aD*(1.0-given) + aC*given), rho)
                    S = wC0 = fitness(q1C, R2*(aC*(1.0-given) + aD*given), rho)
                    if (T<R) and (P<S):
                        # rgb = (0.5-shift, 0.5+shift, 0.5+shift)
                        rgb = (0.5, 0.0, 1.0)
                    else:
                        if (T>R) and (P<S):
                            rgb = (0.0, 1.0, 0.0) # Snowdrift
                        else:
                            rgb = (1.0, 0.0, 0.5)
                            # rgb = (0.5-shift, 0.5+shift, 0.5+shift)
                    ax.plot(aC, aD, color=rgb, marker='o', markerfacecolor=rgb, markersize=3)
                    ax.set_title(title, pad=10.0, fontsize=fstitle)
                    ax.tick_params(axis='x', labelsize=fstick)
                    ax.tick_params(axis='y', labelsize=fstick)
                    ax.set_xlim(0.0, 1.0)
                    ax.set_ylim(0.0, 1.0)
                    ax.set_box_aspect(1)
                    ax.text(0.01, 1.1, letter, fontsize=18, weight='bold')

fig.savefig(filename, transparent=False)
plt.close()

