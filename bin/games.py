#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Parameters

alpha = 0.5
R1 = 2.0
R2 = 2.0
#deathrate = pow(2, -7)
deathrate = 1.0
width = 15
height = 15

n_data = 10
fs = 14

ess = np.linspace(-5, 5, num=11)
ess = pow(2, ess)
rhos = 1.0 - 1.0/ess
givens = np.linspace(1.0, 0.0, num=11)
aDs = np.linspace(0.1, 0.4, 4)
aCs = np.linspace(0.2, 0.5, 4)
colors = ('0.800', '0.600', '0.400', '0.200')
xs = np.linspace(0.0, 1.0, num=n_data)
zerodifs = [0.0] * n_data
repeats = 1.0/(1.0 - pow(1.0 - deathrate, 2))

def fitness(q1, q2, rho):
    if rho == 0.0:
        w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

for aC, aD, color in zip(aCs, aDs, colors):
    fig, axs = plt.subplots(11, 11, figsize=(12, 12))
    q1C = R1*(1.0-aC);
    q1D = R1*(1.0-aD);
    for rowax, given in zip(axs, givens):
        for ax, rho in zip(rowax, rhos):
            wR = fitness(q1C, R2*aC, rho)
            wT = fitness(q1D, R2*(aD*(1-given)+aC*given), rho)
            wS = fitness(q1C, R2*(aC*(1-given)+aD*given), rho)
            wP = fitness(q1D, R2*aD, rho)
            wD = wT*xs + wP*(1-xs)
            wC = wR*xs + wS*(1-xs)
            #difs = wC - wD
            #ax.plot(xs, zerodifs, c='orange')
            #ax.plot(xs, difs, c=color)
            ax.plot(xs, wD, c='orange')
            ax.plot(xs, wC, c='green')
            #ax.set_ylim(-0.1, 0.3)
            ax.set_ylim(0.0, 1.3)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_box_aspect(1)
            #ax.set_xticks(np.linspace(0.0, max_a2, num=3))
            #ax.set_yticks(np.linspace(0.0, 1.5, num=4))
            #ax.xlabel('$\it{a_{2i}}$', fontsize=fs) 
            #ax.ylabel('$\it{w}$', fontsize=fs)

    plt.subplots_adjust(wspace=0, hspace=0)

    plt.savefig(f'games{aC}.png')

    plt.close()
