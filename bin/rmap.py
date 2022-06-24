#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Parameters

alpha = 0.5
rho = 0.0
R1 = 2.0
R2 = 2.0
given = 1.0
deathrate = pow(2, -7)

n_curves = 11
fs = 14

def fitness(a, b):
    if rho == 0.0:
        w = pow(a, alpha)*pow(b, 1.0 - alpha)
    else:
        w = pow(alpha*pow(a, rho) + (1.0 - alpha)*pow(b, rho), 1.0/rho)
    return w

ais = np.linspace(0.0, 1.0, num=n_curves)
repeats = 1.0/(1.0 - pow(1.0 - deathrate, 2))

plt.figure()

for ai in ais:
    y = fitness(R1*(1.0-ai), R2*(ai*(1.0-given) + ais*given))
    plt.plot(ais, y, c='blue')
    y = (fitness(R1*(1.0-ais), R2*(ai*(1.0-given) + ais*given)) + fitness((1.0-ais)*R1, R2*(ais*(1.0-given) + ais*given))*(repeats - 1.0))/repeats
    plt.plot(ais, y, c='red')
    #ax.set_xlim(0.0, max_a2)
    #ax.set_ylim(0.0, 1.5, 1)
    #ax.set_xticks(np.linspace(0.0, max_a2, num=3))
    #ax.set_yticks(np.linspace(0.0, 1.5, num=4))
    plt.xlabel('$\it{a_{2i}}$', fontsize=fs) 
    plt.ylabel('$\it{w}$', fontsize=fs)

# Save

plt.savefig('l.png')

