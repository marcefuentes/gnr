#!/usr/bin/env python

import numpy as np
import pandas as pd

# Parameters

alpha = 0.5
R1 = 2.0
R2 = 2.0
deathrate = 1.0
width = 15
height = 5

n_data = 2
fs = 14

ess = np.linspace(-5, 5, num=11)
ess = pow(2, ess)
rhos = 1.0 - 1.0/ess
givens = np.linspace(1.0, 0.0, num=11)
aC = 0.2
aD = 0.1
q1C = R1*(1.0-aC);
q1D = R1*(1.0-aD);

def fitness(q1, q2, rho):
    if rho == 0.0:
        w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

optimal= []
none = []

for rho in rhos:
    wR = fitness(q1C, R2*aC, rho)
    for given in givens:
        wS = fitness(q1C, R2*(aC*(1.0-given)+aD*given), rho)
        optimal.append([1/(1-rho), given, 1, wR - wS, 0.0])
        none.append([1/(1-rho), given, 1, 0.0, 0.0])

df = pd.DataFrame(optimal, columns=['ES', 'Given', 'Time', 'wmedian', 'wmedianSD'])
df.to_csv('optimal/optimal.csv', index=False)

df = pd.DataFrame(none, columns=['ES', 'Given', 'Time', 'wmedian', 'wmedianSD'])
df.to_csv('none/none.csv', index=False)
