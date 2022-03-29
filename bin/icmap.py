#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# Error message

if len(sys.argv) < 1:
    print('You must enter an argument: ces or q')
    exit(1)

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

# Base equations

MRT = max_a2*R2/(max_a1*R1)

if sys.argv[1] == 'ces':
    Q = R2*pow(MRT*(1.0 - alpha)/alpha,-1.0/(1.0 - rho))/R1
    a1 = max_a2*Q/(1.0 + max_a2*Q/max_a1)   # CES
else:
    a1 = min(max_a1, pow(MRT*c2/(alpha_q*c1), 1/(alpha_q - 1.0))/R1)

q1 = a1*R1
a2 = max_a2 - max_a2*a1/max_a1
q2 = a2*R2

# Fitness functions

def fitness(a, b):
    if sys.argv[1] == 'ces':
        if rho == 0.0:
            w = pow(a, alpha)*pow(b, 1.0 - alpha)
        else:
            w = pow(alpha*pow(a, rho) + (1.0 - alpha)*pow(b, rho), 1.0/rho)
    else:
        w = c1*pow(a, alpha_q) + c2*b
    return w

# x-axis data

x = np.linspace(0.000001, 3.0, num=n)

# Indifference curves

def ic(w):
    if sys.argv[1] == 'ces':
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
    else:
        indiff = (w - c1*pow(x, alpha_q))/c2
    return indiff

# Plot indifference curves and budget line

fs = 16      # Label font size

plt.figure(figsize=(6.5, 6))
ax = plt.gca()
ax.set_xlim(0.0, max_a1*R1*1.5)
ax.set_ylim(0.0, max_a2*R2*1.5)
ax.set_xticks(np.linspace(0.0, max_a1*R1*1.5, num=4))
ax.set_yticks(np.linspace(0.0, max_a2*R2*1.5, num=4))
ax.tick_params(axis='x', labelsize=fs)
ax.tick_params(axis='y', labelsize=fs)
#ax.set_xlabel('$\it{q_1}$', fontsize=fs) 
#ax.set_ylabel('$\it{q_2}$', fontsize=fs)
ax.set_xlabel('Quantity of resource $\it{A}$', fontsize=fs) 
ax.set_ylabel('Quantity of resource $\it{B}$', fontsize=fs) 

for w in np.linspace(0.4, 1.6, num=n_ic):
    plt.plot(x, ic(w), c='#dbdbdb')

max_w = fitness(q1, q2)
plt.plot(x, ic(max_w), c='#7e7e7e')

budget = max_a2*R2 - MRT*x
plt.plot(x, budget, c='black')

# Save

plt.savefig('ICmap' + sys.argv[1] + '.png')

# q1

x = np.linspace(0.0, max_a2, num=n)

# Plot fitness landscape

plt.figure()
ax = plt.gca()
ax.set_xlim(0.0, max_a2)
ax.set_ylim(0.0, 1.5, 1)
ax.set_xticks(np.linspace(0.0, max_a2, num=3))
ax.set_yticks(np.linspace(0.0, 1.5, num=4))
ax.set_xlabel('$\it{a_2}$', fontsize=fs) 
ax.set_ylabel('$\it{w}$', fontsize=fs)
plt.plot(x, fitness((max_a2 - x)*max_a1*R1/max_a2, x*R2), c='black')

# Save

plt.savefig('landscape' + sys.argv[1] + '.png')

