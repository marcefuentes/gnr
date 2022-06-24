#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Parameters

alpha = 0.5
ES = pow(2, -5)
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

rho = 1.0 - 1.0/ES
ais = np.linspace(0.0, 1.0, num=n_curves)
repeats = 1.0/(1.0 - pow(1.0 - deathrate, 2))

plt.figure()

for ai in ais:
    y0 = fitness(R1*(1.0-ais), R2*(ais*(1.0-given) + ai*given))
    plt.plot(ais, y0, c='yellow', linewidth=ai*2.0)
    yr = fitness(R1*(1.0-ai), R2*(ai*(1.0-given) + ai*given))
    y = (y0 + yr*(repeats - 1.0))/repeats
    plt.plot(ais, y, c='red', linewidth=ai*2.0)
    y = ((y0 + yr)/2.0)
    #plt.plot(ais, y, c='orange', linewidth=ai*2.0)
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 2.0)
    #ax.set_xticks(np.linspace(0.0, max_a2, num=3))
    #ax.set_yticks(np.linspace(0.0, 1.5, num=4))
    plt.xlabel('$\it{a_{2i}}$', fontsize=fs) 
    plt.ylabel('$\it{w}$', fontsize=fs)

# Save

plt.savefig('l.png')

