#! /usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.colors import ListedColormap
import numpy as np
import time

start_time = time.perf_counter ()

alpha = 0.75
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0
npoints = 128

num = 21    # Number of subplot rows and columns
every = int(num/2)
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

def fitness(x, y, given, rho):
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - given) + x*R2*given
    if rho == 0.0:
        w = q1*q2
        if w > 0.0: w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
    elif rho < 0.0:
        w = q1*q2
        if w > 0.0: w = alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho)
        if w > 0.0: w = pow(w, 1.0/rho)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

R = R2/R1
b = a2max/a1max
givens = np.linspace(maxgiven, mingiven, num=num)
givens[0] = 0.9999999
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
RR, GG = np.meshgrid(rhos, givens)

fslabel = 26 # Label font size
fstick = 18 # Tick font size

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(12, 6), constrained_layout=False) 
fig.supylabel("Partner's share of $\it{A}$", x=0.05, y=0.520, fontsize=fslabel)
fig.supxlabel("Substitutability of $\it{A}$", x=0.525, fontsize=fslabel)

outer_grid = fig.add_gridspec(1, 2, left=0.15, right=0.9, top=0.86, bottom=0.176)

# Continuous

grid = outer_grid[0, 0].subgridspec(num, num, wspace=0, hspace=0)
axs = grid.subplots()

a2s = np.linspace(0.0, a2max, num=npoints)
TT = b*R*(1.0 - GG)
Q = R*pow(TT*(1.0 - alpha)/alpha, 1.0/(RR - 1.0))
a2eqss = a2max/(1.0 + Q*b)

for row, given, a2eqs in zip(axs, givens, a2eqss):
    for ax, rho, a2eq in zip(row, rhos, a2eqs):
        w = []
        for a2 in a2s:
            w.append(fitness(a2eq, a2, given, rho))
        weq = fitness(a2eq, a2eq, given, rho)
        ax.plot(a2s, w, linewidth=2, c=cm.magma(weq))
        ax.set(xticks=[], yticks=[], xlim=(0.0, a2max), ylim=(0.0, 2.0))
        ax.set_box_aspect(1)
axs[0, 0].set_title('a', fontsize=fslabel, weight='bold')
for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
    ax.set_xlabel(round(log_es), fontsize=fstick)
for ax, given in zip(axs[::every, 0], givens[::every]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

# Discrete

grid = outer_grid[0, 1].subgridspec(num, num, wspace=0, hspace=0)
axs = grid.subplots()

a2 = np.linspace(0.0, a2max, num=3)

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        T = fitness(a2[1], a2[0], given, rho)
        R = fitness(a2[1], a2[1], given, rho)
        P = fitness(a2[0], a2[0], given, rho)
        S = fitness(a2[0], a2[1], given, rho)
        x = 0.0
        if (T < R) & (P < S): x = 1.0
        elif (T > R) & (P < S) & (R - S - T + P != 0.0): x = (P - S)/(R - S - T + P)
        weq = (T + S)*x*(1.0 - x) + R*x*x + P*(1.0 - x)*(1.0 - x)
        w = []
        for a in a2:
            w.append(fitness(a2[0], a, given, rho)*(1.0 - x) + fitness(a2[1], a, given, rho)*x)   
        ax.plot(a2, w, linewidth=2, c=cm.magma(weq))
        ax.set(xticks=[], yticks=[], xlim=(0.0, a2max), ylim=(0.0, 2.0))
        ax.set_box_aspect(1)
axs[0, 0].set_title('b', fontsize=fslabel, weight='bold')
for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
    ax.set_xlabel(round(log_es), fontsize=fstick)

plt.savefig('map.png', dpi=100)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
