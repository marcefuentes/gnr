#!/usr/bin/env python

import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

width = 4
height = 4
fslabel = 15
fstick = 10    

# Parameters

alpha = 0.5
R1 = 2.0
R2 = 2.0

ess = np.linspace(-5, 5, num=11)
ess = pow(2, ess)
rhos = 1.0 - 1.0/ess
givens = np.linspace(1.0, 0.0, num=11)
aC = 0.3
aD = 0.2
q1C = R1*(1.0-aC);
q1D = R1*(1.0-aD);

def fitness(q1, q2, rho):
    if rho == 0.0:
        w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

x = []
y = []
size = []
color = []

for rho in rhos:
    R = wC1 = fitness(q1C, R2*aC, rho)
    P = wD0 = fitness(q1D, R2*aD, rho)
    for given in givens:
        T = wD1 = fitness(q1D, R2*(aD*(1.0-given) + aC*given), rho)
        S = wC0 = fitness(q1C, R2*(aC*(1.0-given) + aD*given), rho)
        diff = (P-S)*(R-P)
        if diff > 0:
            color.append('#05a630') 
        else:
            color.append('#c9beec')
        size.append(abs(diff))
        x.append([1/(1-rho)])
        y.append(given)

fig, ax = plt.subplots(figsize=(width,height))
ax.scatter(x=x, y=y, s=16000.0*np.array(size), color=color, ec=color)
ax.set_xlabel('Substitutability of resource $\it{A}$', fontsize=fslabel)
ax.set_ylabel('Partner\'s share of resource $\it{A}$', fontsize=fslabel)
ax.tick_params(axis='x', labelsize=fstick)
ax.tick_params(axis='y', labelsize=fstick)
ax.set_xscale('log', base=2)
ax.set_box_aspect(1)

fig.savefig('games.png', bbox_inches='tight', transparent=False)
plt.close()

