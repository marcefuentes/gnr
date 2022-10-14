#! /usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.colors import ListedColormap
import numpy as np
import time

start_time = time.perf_counter ()

alpha = 0.5
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0
npoints = 128
npoints_ic = 128

num = 5     # Number of subplot rows and columns
n_ic = 5    # Number of indifference curves
every = int(num/2)
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

def fitness(x, y, given, rho):
    q1 = (1.0 - y)*R1
    q2 = y*R2*(1.0 - given) + x*R2*given
    if rho == 0.0:
        w = q1*q2
        w[w>0.0] = pow(q1[w>0.0], alpha)*pow(q2[w>0.0], 1.0 - alpha)
    elif rho < 0.0:
        w = q1*q2
        w[w>0.0] = alpha*pow(q1[w>0.0], rho) + (1.0 - alpha)*pow(q2[w>0.0], rho)
        w[w>0.0] = pow(w[w>0.0], 1.0/rho)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

def icces(q, w, rho):
    if rho == 0.0:
        q2 = np.piecewise(q, [q == 0.0, q > 0.0], [1000.0, lambda i: pow(w/pow(i, alpha), 1.0/(1.0 - alpha))])
    elif rho < 0.0:
        q2 = np.piecewise(q, [q == 0.0, q > 0.0], [1000.0, lambda i: np.piecewise(i, [pow(w, rho) <= alpha*pow(i, rho), pow(w, rho) > alpha*pow(i, rho)], [1000.0, lambda j: pow((pow(w, rho) - alpha*pow(j, rho))/(1.0 - alpha), 1.0/rho)])])
    else:
        q2 = np.piecewise(q, [pow(w, rho) <= alpha*pow(q, rho), pow(w, rho) > alpha*pow(q, rho)], [-0.1, lambda i: pow((pow(w, rho) - alpha*pow(i, rho))/(1.0 - alpha), 1.0/rho)])
    return q2

R = R2/R1
b = a2max/a1max
q1b = np.array([0.0, a1max*R1/2.0, a1max*R1])
q1 = np.linspace(0.0, a1max*R1, num=npoints_ic)
ws = np.linspace(2.0/(n_ic + 1), 2.0*n_ic/(n_ic + 1), num=n_ic)
givens = np.linspace(maxgiven, mingiven, num=num)
givens[0] = 0.9999999
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)

fslabel = 26 # Label font size
fstick = 18 # Tick font size

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(12, 6), constrained_layout=False) 
fig.supylabel("Partner's share of $\it{A}$", y=0.520, fontsize=fslabel)
fig.supxlabel("Substitutability of $\it{A}$", x=0.525, fontsize=fslabel)

outer_grid = fig.add_gridspec(1, 2, left=0.15, right=0.9, top=0.86, bottom=0.176)

# Indifference curves and continuous budget line

grid = outer_grid[0, 0].subgridspec(num, num, wspace=0, hspace=0)
axs = grid.subplots()

Ts = b*R*(1.0 - givens)
RR, TT = np.meshgrid(rhos, Ts)
Q = R*pow(TT*(1.0 - alpha)/alpha, 1.0/(RR - 1.0))
a2eqs = a2max/(1.0 + Q*b)

for row, given, g in zip(axs, givens, a2eqs):
    for ax, rho, a2 in zip(row, rhos, g):
        for w in ws:
            ax.plot(q1, icces(q1, w, rho), c='0.950')
        T = b*R*(1.0 - given)
        budget = a2max*R2*(1.0 - given) + a2*R2*given - T*q1b
        ax.plot(q1b, budget, c='green')
        a = np.array([a2])
        weq = fitness(a, a, given, rho)
        ax.plot(q1, icces(q1, weq[0], rho), c=cm.magma(weq))
        ax.set(xticks=[], yticks=[], xlim=(0, R1*a1max), ylim=(0, R2*a2max))
        ax.set_box_aspect(1)
axs[0, 0].set_title('a', fontsize=fslabel, weight='bold')
for ax, given in zip(axs[::every, 0], givens[::every]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)
for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
    ax.set_xlabel(round(log_es), fontsize=fstick)

# Indifference curves and discrete budget line

grid = outer_grid[0, 1].subgridspec(num, num, wspace=0, hspace=0)
axs = grid.subplots()

a2 = np.array([0.5, 0.0])
X, Y = np.meshgrid(a2, a2)

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        for w in ws:
            ax.plot(q1, icces(q1, w, rho), c='0.950')
        Z = fitness(X, Y, given, rho)
        T = Z[1, 0]
        R = Z[0, 0]
        P = Z[1, 1]
        S = Z[0, 1]
        if (T < R) and (P < S):
            xeq = 1.0
        elif (T > R) and (P >= S):
            xeq = 0.0
        else:
            xeq = (P - S)/(R - S - T + P) 
        a2eq = xeq*a2[0]
        T = b*R*(1.0 - given)
        budget = a2max*R2*(1.0 - given) + a2eq*R2*given - T*q1b
        ax.plot(q1b, budget, c='green', marker='o', markersize=5, linestyle='dashed')
        weq = (T + S)*xeq*(1.0 - xeq) + R*xeq*xeq + P*(1.0 - xeq)*(1.0 - xeq)
        ax.plot(q1, icces(q1, weq, rho), c=cm.magma(weq))
        ax.set(xticks=[], yticks=[], xlim=(0, R1*a1max), ylim=(0, R2*a2max))
        ax.set_box_aspect(1)
axs[0, 0].set_title('b', fontsize=fslabel, weight='bold')
for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
    ax.set_xlabel(round(log_es), fontsize=fstick)

plt.savefig('icurves.png', dpi=100)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
