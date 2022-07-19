#!/usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = f'gamesTR{sys.argv[1]}.png'

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

width = 11.0
height = 6.0
fslabel = 26
fstitle= 24
fstick = 16    

# Parameters

alpha = 0.5
R1 = 2.0
R2 = 2.0

ess = np.linspace(-5, 5, num=11) if sys.argv[1] == 'ces' else np.linspace(-10, 0, num=11)
ess = pow(2, ess)
givens = np.linspace(1.0, 0.0, num=11)
aCs = [0.3, 0.5]
aDs = [0.2, 0.4]

def fitness(q1, q2, rho):
    if sys.argv[1] == 'ces':
        if rho == 0.0:
            w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
        else:
            w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    else:
        w = 4.0*pow(q1, rho)/9.0 + 4.0*q2/9.0
    return w

def dif_color(dif):
    color = (red+4.0*dif, green+4.0*dif, blue) if dif <= 0.0 else (red-4.0*dif, green, blue-4.0*dif)
    return color

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
            if (T<R) and (P<S):
                rgb_edge = (0.97, 0.97, 0.97) # No dilemma
                rgb = (0.97, 0.97, 0.97)
            elif (T>R) and (P<S):
                rgb_edge = (0.8, 1.0, 0.8) # Snowdrift
                rgb = (1.0-8*(T-R), 1.0-4*(T-R), 1.0-8*(T-R))
            else:
                rgb_edge = (0.8, 0.8, 1.0) # Prisoner's dilemma
                rgb = (1.0-8*(T-R), 1.0-4*(T-R), 1.0-8*(T-R))
            #rgb = (0.7-T+R, 1.0, 0.7+P-S)
            color.append(rgb)
            edgecolor.append(rgb)
            size.append(1200.0*abs(R-P))

fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(width,height), constrained_layout=False, squeeze=False)
axs = axs.flatten()
fig.supxlabel('Substitutability of $\it{A}$', fontsize=fslabel)
fig.supylabel('Partner\'s share of $\it{A}$', fontsize=fslabel)

titles = ['$\it{T}$ - $\it{R}$ and $\it{R}$ - $\it{P}$\nfor $\it{a}$ = {0.2, 0.3}', '$\it{T}$ - $\it{R}$ and $\it{R}$ - $\it{P}$\nfor $\it{a}$ = {0.4, 0.5}']
for size, color, edgecolor, title, ax in zip(sizes, colors, edgecolors, titles, axs): 
    ax.scatter(x=x, y=y, s=size, color=color, ec=edgecolor)
    ax.set_title(title, pad=10.0, fontsize=fstitle)
    ax.tick_params(axis='x', labelsize=fstick)
    ax.tick_params(axis='y', labelsize=fstick)
    ax.set_xscale('log', base=2)
    ax.set_box_aspect(1)

fig.savefig(filename, transparent=False)
plt.close()

