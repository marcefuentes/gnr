#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys

filename = f'maxw{sys.argv[1]}.png'

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

height = 6.0
fslabel = 48
fsnumbers= 24

# Parameters

alpha = 0.5
R1 = 2.0
R2 = 2.0
npoints = 32

log_ess = np.linspace(-5, 5, num=11) if sys.argv[1] == 'ces' else np.linspace(-10, 0, num=11)
ess = pow(2, log_ess)
givens = np.linspace(1.0, 0.0, num=11)
x = np.linspace(1/(npoints + 1), 1.0 - 1/(npoints + 1), num = npoints)
y = np.linspace(1.0 - 1/(npoints + 1), 1/(npoints + 1), num = npoints)
X, Y = np.meshgrid(x, y)

def fitness(x, y, given, rho):
    q1 = (1-y)*R1
    q2 = y*R2*(1-given) + x*R2*given
    if sys.argv[1] == 'q':
        w = 4.0*pow(q1, rho)/9.0 + 4.0*q2/9.0
    else:
        if rho == 0.0:
            w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
        else:
            w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

fig, axs = plt.subplots(nrows=11, ncols=11, figsize=(14.0, 14.0), constrained_layout=False)
fig.supxlabel('Substitutability of $\it{A}$', fontsize=fslabel, x=0.513)
fig.supylabel('Partner\'s share of $\it{A}$', fontsize=fslabel)

for row, given in zip(axs, givens):
    for ax, es, log_es in zip(row, ess, log_ess):
        rho = 1.0 - 1.0/es if sys.argv[1] == 'ces' else es
        Z = fitness(X, Y, given, rho)
        nZ = Z
        nZ[np.argmax(Z, axis=0), np.arange(Z.shape[1])] = 2
        ax.imshow(nZ, vmin=0, vmax=2)
        ax.set(xticks=[], yticks=[], xlim=(0, npoints-2), ylim=(npoints-2, 0))
        #ax.set_box_aspect(1)
        #ax.text(0.01, 1.1, letter, fontsize=18, weight='bold')
        if given == 0.0:
            x = '$2^{{{}}}$'.format(round(log_es))
            ax.set_xlabel(x, fontsize=fsnumbers)
        if log_es == -5:
            ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fsnumbers)

fig.savefig(filename, transparent=False)
plt.close()

