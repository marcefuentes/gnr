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

# Parameters

alpha = 0.5
R1 = 2.0
R2 = 2.0
mutationsize = 0.0078125

#log_ess = np.linspace(-5, 5, num=11) if sys.argv[1] == 'ces' else np.linspace(-10, 0, num=11)
log_ess = np.linspace(0, 0, num=1) if sys.argv[1] == 'ces' else np.linspace(-10, 0, num=11)
ess = pow(2, log_ess)
#givens = np.linspace(1.0, 0.0, num=11)
givens = [ 1.0 ]
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
            aC = mutationsize*2
            while aC < 1.0:
                q1C = R1*(1.0-aC);
                R = wC1 = fitness(q1C, R2*aC, rho)
                aD = mutationsize
                while aD < aC:
                    q1D = R1*(1.0-aD);
                    P = wD0 = fitness(q1D, R2*aD, rho)
                    T = wD1 = fitness(q1D, R2*(aD*(1.0-given) + aC*given), rho)
                    S = wC0 = fitness(q1C, R2*(aC*(1.0-given) + aD*given), rho)
                    if R>P:
                        red = 1.0/(1.0 + math.exp((T-R-0.4)*10.0))
                        blue = 1.0/(1.0 + math.exp((P-S-0.4)*10.0))
                    else:
                        red = 1.0/(1.0 + math.exp(-(T-R-0.4)*10.0))
                        blue = 1.0/(1.0 + math.exp(-(P-S-0.4)*10.0))
                    rgb = (red, 0.9, blue) # No dilemma
                    if (T<R) and (P<S):
                        #rgb = (0.5, 0.0, 1.0)
                        rgb = (red, 0.9, blue) # No dilemma
                    else:
                        if (T>R) and (P<S):
                            #rgb = (0.0, 0.1, 0.0)
                            rgb = (red, 0.9, blue) # Snowdrift
                        else:
                            #rgb = (1.0, 0.0, 0.5)
                            rgb = (red, 0.9, blue) # Prisoner's dilemma
                    ax.plot(aC, aD, color=rgb, marker='o', markerfacecolor=rgb, markersize=math.fabs((R-P)*1.0))
                    #ax.set_title(title, pad=10.0, fontsize=fstitle)
                    #ax.tick_params(axis='x', labelsize=fstick)
                    #ax.tick_params(axis='y', labelsize=fstick)
                    #ax.set_xlim(0.0, 1.0)
                    #ax.set_ylim(0.0, 1.0)
                    #ax.set_box_aspect(1)
                    #ax.text(0.01, 1.1, letter, fontsize=18, weight='bold')
                    aD += mutationsize
                aC += mutationsize

fig.savefig(filename, transparent=False)
plt.close()

