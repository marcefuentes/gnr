#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

fs = 16      # Label font size

# Parameters

alpha = 0.5
rho = 0.0
R1 = 2.0
R2 = 2.0
max_a1 = 1.0
max_a2 = 1.0

alpha_q = 0.125
c1 = 4/9
c2 = 4/9

n = 100     # Number of x-axis points
n_ic = 7    # Number of indifference curves

# Indifference curves

def icces(w):
    if rho == 0.0:
        indiff = pow(w/pow(x, alpha), 1.0/(1.0 - alpha))
    else:
        indiff = [0.0]*n
        for i in range(n):
            if pow(w, rho) <= alpha*pow(x[i], rho):
                if rho < 0.0:
                    indiff[i] = 1000.0
                else:
                    indiff[i] = -0.1
            else:
                indiff[i] = pow((pow(w, rho) - alpha*pow(x[i], rho))/(1.0 - alpha), 1.0/rho)
    return indiff

def icquasilinear(w):
    indiff = (w - c1*pow(x, alpha_q))/c2
    return indiff

# Fitness landscapes

def wces(a, b):
    if rho == 0.0:
        w = pow(a, alpha)*pow(b, 1.0 - alpha)
    else:
        w = pow(alpha*pow(a, rho) + (1.0 - alpha)*pow(b, rho), 1.0/rho)
    return w

def wquasilinear(a, b):
    w = c1*pow(a, alpha_q) + c2*b
    return w

MRT = max_a2*R2/(max_a1*R1)

# For ces

Q = R2*pow(MRT*(1.0 - alpha)/alpha,-1.0/(1.0 - rho))/R1
a1 = max_a2*Q/(1.0 + max_a2*Q/max_a1)   # CES
q1 = a1*R1
a2 = max_a2 - max_a2*a1/max_a1
q2 = a2*R2
cesdict = {'wfunction': wces, 'icfunction': icces, 'q1': q1, 'q2': q2}

# For quasilinear

a1 = min(max_a1, pow(MRT*c2/(alpha_q*c1), 1/(alpha_q - 1.0))/R1)
q1 = a1*R1
a2 = max_a2 - max_a2*a1/max_a1
q2 = a2*R2
qdict = {'wfunction':  wquasilinear, 'icfunction': icquasilinear, 'q1': q1, 'q2': q2}

ufs = [cesdict, qdict]

# x-axis data

x = np.linspace(0.000001, 3.0, num=n)

# Plot indifference curves and budget line

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 7), sharey=True)
fig.supxlabel('Quantity of resource $\it{A}$ consumed', fontsize=fs+4) 
fig.supylabel('Quantity of resource $\it{B}$ consumed', x=0.02, fontsize=fs+4, ha='center') 

for ax, uf in zip(axs, ufs):
    ax.set_xlim(0.0, max_a1*R1*1.5)
    ax.set_ylim(0.0, max_a2*R2*1.5)
    ax.set_xticks(np.linspace(0.0, max_a1*R1*1.5, num=4))
    ax.set_yticks(np.linspace(0.0, max_a2*R2*1.5, num=4))
    ax.tick_params(axis='x', labelsize=fs)
    ax.tick_params(axis='y', labelsize=fs)

    for w in np.linspace(0.4, 1.6, num=n_ic):
        ax.plot(x, uf['icfunction'](w), c='#dbdbdb')

    max_w = uf['wfunction'](uf['q1'], uf['q2'])
    ax.plot(x, uf['icfunction'](max_w), c='#7e7e7e')

    budget = max_a2*R2 - MRT*x
    ax.plot(x, budget, c='black')
    ax.set_box_aspect(1)

plt.savefig('ICmap.png', transparent=False)

# q1

x = np.linspace(0.0, max_a2, num=n)

# Plot fitness landscape

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 7), sharey=True)
fig.supxlabel('Effort allocated to capturing resource $\it{B}$', fontsize=fs+4) 
fig.supylabel('Fitness', fontsize=fs+4, x=0.02, ha='center')

for ax, uf in zip(axs, ufs):
    ax.set_xlim(0.0, max_a2)
    ax.set_ylim(0.0, 1.5, 1)
    ax.set_xticks(np.linspace(0.0, max_a2, num=3))
    ax.set_yticks(np.linspace(0.0, 1.5, num=4))
    ax.tick_params(axis='x', labelsize=fs)
    ax.tick_params(axis='y', labelsize=fs)
    ax.plot(x, uf['wfunction']((max_a2 - x)*max_a1*R1/max_a2, x*R2), c='black')
    ax.set_box_aspect(1)

plt.savefig('landscape.png')

