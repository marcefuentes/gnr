#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Parameters

alpha = 0.5
R1 = 2.0
R2 = 2.0
deathrate = pow(2, -7)
width = 15
height = 15

n_curves = 11
fs = 14

ess = np.linspace(-5, 5, num=11)
ess = pow(2, ess)
rhos = 1.0 - 1.0/ess
givens = np.linspace(1.0, 0.0, num=11)
ais = np.linspace(0.0, 1.0, num=n_curves)
repeats = 1.0/(1.0 - pow(1.0 - deathrate, 2))

def fitness(q1, q2, rho):
    if rho == 0.0:
        w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

fig, axs = plt.subplots(11, 11, figsize=(12, 12))

for rowax, given in zip(axs, givens):
    for ax, rho in zip(rowax, rhos):
        for ai in ais:
            y05 = fitness(R1*(1.0-ai), R2*(ai*(1.0-given) + 0.5*given), rho)
            ax.scatter(ai, y0, c='green', s=3)
            yr = fitness(R1*(1.0-0.5), R2*(0.5*(1.0-given) + 0.5*given), rho)
            y = (y0 + yr*(repeats - 1.0))/repeats
            #ax.scatter(ai, y1, c='yellow', s=3)
            y = ((y0 + yr)/2.0)
            ax.scatter(ai, y, c='orange', s=3)
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, 2.0)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_box_aspect(1)
            #ax.set_xticks(np.linspace(0.0, max_a2, num=3))
            #ax.set_yticks(np.linspace(0.0, 1.5, num=4))
            #ax.xlabel('$\it{a_{2i}}$', fontsize=fs) 
            #ax.ylabel('$\it{w}$', fontsize=fs)

plt.subplots_adjust(wspace=0, hspace=0)

plt.savefig('l.png')
plt.close()
