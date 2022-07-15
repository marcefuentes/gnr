#!/usr/bin/env python

import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

width = 4
height = 4
fslabel = 15
fstick = 10    
red = 0.97
green = 0.97
blue = 0.97

# Parameters

alpha = 0.5
R1 = 2.0
R2 = 2.0

ess = np.linspace(-5, 5, num=11) if sys.argv[1] == 'ces' else np.linspace(-10, 0, num=11)
ess = pow(2, ess)
givens = np.linspace(1.0, 0.0, num=11)
aC = 0.5
aD = 0.4
q1C = R1*(1.0-aC);
q1D = R1*(1.0-aD);

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
size = []
color = []

for es in ess:
    rho = 1.0 - 1.0/es if sys.argv[1] == 'ces' else es
    R = wC1 = fitness(q1C, R2*aC, rho)
    P = wD0 = fitness(q1D, R2*aD, rho)
    for given in givens:
        T = wD1 = fitness(q1D, R2*(aD*(1.0-given) + aC*given), rho)
        S = wC0 = fitness(q1C, R2*(aC*(1.0-given) + aD*given), rho)
        if (T<R) and (P<S):
            rgb = (red-0.1, green-0.4-R+T, blue-0.4+P-S)    # No dilemma
        elif (T>R) and (P<S):
            rgb = (red-0.4-T+R, green-0.1, blue-0.4+P-S)    # Snowdrift
        else:
            rgb = (red-0.4-T+R, green-0.4-P+S, blue-0.1)    # Prisoner's dilemma
        color.append(rgb)
        size.append(abs(R-P))
        x.append(es)
        y.append(given)

fig, ax = plt.subplots(figsize=(width,height))
ax.scatter(x=x, y=y, s=600.0*np.array(size), color=color, ec=color)
ax.set_xlabel('Substitutability of resource $\it{A}$', fontsize=fslabel)
ax.set_ylabel('Partner\'s share of resource $\it{A}$', fontsize=fslabel)
ax.tick_params(axis='x', labelsize=fstick)
ax.tick_params(axis='y', labelsize=fstick)
ax.set_xscale('log', base=2)
ax.set_box_aspect(1)

fig.savefig(f'games{sys.argv[1]}.png', bbox_inches='tight', transparent=False)
plt.close()

