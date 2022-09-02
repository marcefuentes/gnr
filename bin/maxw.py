#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = f'maxw{sys.argv[1]}.png'

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

height = 6.0
fslabel = 40
fstitle= 30

# Parameters

alpha = 0.5
R1 = 2.0
R2 = 2.0
mutationsize = 0.0078125

log_ess = np.linspace(-5, 5, num=11) if sys.argv[1] == 'ces' else np.linspace(-10, 0, num=11)
ess = pow(2, log_ess)
givens = np.linspace(1.0, 0.0, num=11)
aXs = aYs = np.linspace(mutationsize, 1.0 - mutationsize)

def fitness(q1, q2, rho):
    if sys.argv[1] == 'q':
        w = 4.0*pow(q1, rho)/9.0 + 4.0*q2/9.0
    else:
        if rho == 0.0:
            w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
        else:
            w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

fig, axs = plt.subplots(nrows=11, ncols=11, sharey=True, figsize=(12.0, 12.0), constrained_layout=False, squeeze=False)
fig.supxlabel('Substitutability of $\it{A}$', fontsize=fslabel)
fig.supylabel('Partner\'s share of $\it{A}$', fontsize=fslabel)

for row, given in zip(axs, givens):
    for ax, es, log_es in zip(row, ess, log_ess):
        rho = 1.0 - 1.0/es if sys.argv[1] == 'ces' else es
        for aX in aXs:
            q1X = R1*(1.0-aX)
            for aY in aYs:
                q1Y = R1*(1.0-aY)
                T = fitness(q1Y, R2*(aY*(1.0-given) + aX*given), rho)
                rgb = (T/2.0, T/2.0, T/2.0) 
                ax.plot(aX, aY, color=rgb, marker='s', markerfacecolor=rgb, markersize=3.0)
                ax.set(xticks=[], yticks=[], xlim=[0.0, 1.0], ylim=[0.0, 1.0])
                #ax.set_box_aspect(1)
                #ax.text(0.01, 1.1, letter, fontsize=18, weight='bold')
        if given == 0.0:
            x = '$2^{{{}}}$'.format(log_es)
            ax.set_xlabel(x, fontsize=fslabel)
        if log_es == -5:
            print(given)
            ax.set_ylabel(given, rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fslabel)

fig.savefig(filename, transparent=False)
plt.close()

