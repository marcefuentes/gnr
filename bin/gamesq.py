#!/usr/bin/env python

import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

fs = 16     # Label font size

# Parameters

alpha = 0.5
R1 = 2.0
R2 = 2.0

ess = np.linspace(-10, 0, num=11)
ess = pow(2, ess)
rhos = ess
#rhos = 1.0 - 1.0/ess
givens = np.linspace(1.0, 0.0, num=11)
aC = 0.5
aD = 0.2
q1C = R1*(1.0-aC);
q1D = R1*(1.0-aD);

def fitness(q1, q2, rho):
    if rho == 0.0:
        w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

def fitnessq(q1, q2, rho):
    w = 4.0*pow(q1, rho)/9.0 + 4.0*q2/9.0
    return w

x = []
y = []
size0 = []
size1 = []
color0 = []
color1 = []

for rho in rhos:
    R = fitnessq(q1C, R2*aC, rho)
    P = fitnessq(q1D, R2*aD, rho)
    for given in givens:
        #x.append([1/(1-rho)])
        x.append(rho)
        y.append(given)
        T = fitnessq(q1D, R2*(aD*(1.0-given) + aC*given), rho)
        S = fitnessq(q1C, R2*(aC*(1.0-given) + aD*given), rho)
        wC0 = S
        wD0 = P
        wC1 = R
        wD1 = T
        diff0 = wC0 - wD0
        if diff0 > 0:
            color0.append('#59ff00') 
        else:
            color0.append('#ff8300')
        size0.append(abs(diff0))
        diff1 = wC1 - wD1
        if diff1 > 0:
            color1.append('#2d8000') 
        else:
            color1.append('#ff0000')
        size1.append(abs(diff1))

fig, ax = plt.subplots(figsize=(7,7))
ax.scatter(x=x, y=y, s=1200.0*np.array(size0), color='none', ec=color0, linewidth=3, alpha=0.5)
ax.scatter(x=x, y=y, s=1200.0*np.array(size1), color='none', ec=color1, linewidth=3, alpha=0.5)
ax.set_xlabel('Substitutability of resource $\it{A}$', fontsize=fs+6)
ax.set_ylabel('Partner\'s share of resource $\it{A}$', fontsize=fs+6)
ax.tick_params(axis='x', labelsize=fs)
ax.tick_params(axis='y', labelsize=fs)
ax.set_xscale('log', base=2)
ax.set_box_aspect(1)

fig.savefig('gamesq.png', transparent=False)
plt.close()

