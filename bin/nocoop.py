#! /usr/bin/env python

from math import log
import matplotlib.pyplot as plt
import numpy as np
import time

start_time = time.perf_counter ()

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

letters = [['a', 'b', 'c', 'd', 'e'],
            ['f', 'g', 'h', 'i', 'j'],
            ['k', 'l', 'm', 'n', 'o'],
            ['p', 'q', 'r', 's', 't'],
            ['u', 'v', 'w', 'x', 'y']]

traits = ['a2Seenmedian', 'helpmedian', 'wmedian']
traitlabels = ['Effort to get $\it{A}$', 'Help', 'Fitness']
traitvmaxs = [0.5, 1.0, 1.0]

alpha = 0.5
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

num = 2000
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

def a2eq(X, Y):
    T = b*R*(1.0 - Y)
    Q = R*pow(T*(1.0 - alpha)/alpha, 1.0/(X - 1.0))
    a2 = a2max/(1.0 + Q*b)
    return a2

def fitness(Z, X, Y):
    q1 = (1.0 - Z)*R1
    q2 = Z*R2*(1.0 - Y) + Z*R2*Y
    w = np.where(X == 0.0, pow(q1, alpha)*pow(q2, 1.0 - alpha), pow(alpha*pow(q1, X) + (1.0 - alpha)*pow(q2, X), 1.0/X)) 
    return w

R = R2/R1
b = a2max/a1max
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
givens = np.linspace(maxgiven, mingiven, num=num)
givens[0] = 0.99999
Xrhos, Ygivens = np.meshgrid(rhos, givens)
Z0 = a2eq(Xrhos, Ygivens)
Zs = [Z0, Z0*R2*Ygivens, fitness(Z0, Xrhos, Ygivens)]
Zmaxs = [0.5, 1.0, 1.0]

fslabel=24 # Label font size
fstick=16 # Tick font size

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(18, 6)) # constrained_layout=False, squeeze=False

fig.supxlabel('Substitutability of $\it{A}$', x=0.513, y=0.00, fontsize=fslabel*1.25)
fig.supylabel('Partner\'s share of $\it{A}$', x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

axs.flatten()
extent = 0, num, 0, num
for ax, Z, Zmax, letter, traitlabel in zip(axs, Zs, Zmaxs, letters[0], traitlabels):
    ax.imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=Zmax)
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)
    ax.set_xticks([0, num/2.0, num])
    ax.set_yticks([0, num/2.0, num])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if traitlabel == traitlabels[0]:
        ax.set_yticklabels([round(mingiven, 1), round((mingiven + maxgiven)/2, 1), round(maxgiven, 1)], fontsize=fstick) 
    ax.text(0, 2070, letter, fontsize=fslabel, weight='bold')
    ax.set_xticklabels([round(minlog_es, 1), round((minlog_es + maxlog_es)/2, 1), round(maxlog_es, 1)], fontsize=fstick)

plt.savefig('nocoop.png', transparent=False)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
