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
npoints_ic = 128

num = 5    # Number of subplot rows and columns
#markersize = 10
n_ic = 5    # Number of indifference curves
every = int(num/2)
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

def fitness(x, y):
    q1 = (1.0 - y)*R1
    q2 = y*R2*(1.0 - GG) + x*R2*GG
    w = q1*q2
    mask = (w > 0.0) & (RR == 0.0)
    w[mask] = pow(q1[mask], alpha)*pow(q2[mask], 1.0 - alpha)
    mask = (w > 0.0) & (RR < 0.0)
    w[mask] = alpha*pow(q1[mask], RR[mask]) + (1.0 - alpha)*pow(q2[mask], RR[mask])
    mask = (w > 0.0) & (RR < 0.0)
    w[mask] = pow(w[mask], 1.0/RR[mask])
    mask = (RR > 0.0)
    w[mask] = pow(alpha*pow(q1[mask], RR[mask]) + (1.0 - alpha)*pow(q2[mask], RR[mask]), 1.0/RR[mask])
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
a1_budget = np.linspace(0.0, a1max, num=3)
q2_budget = (a2max - b*a1_budget)*R2
q1_budget = a1_budget*R1
q1_ic = np.linspace(0.0, a1max*R1, num=npoints_ic)
ws = np.linspace(2.0/(n_ic + 1), 2.0*n_ic/(n_ic + 1), num=n_ic)
givens = np.linspace(maxgiven, mingiven, num=num)
givens[0] = 0.9999999
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
RR, GG = np.meshgrid(rhos, givens)
icss = []
for rho in rhos:
    ics = []
    for w in ws:
        ics.append(icces(q1_ic, w, rho))
    icss.append(ics)

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

TT = b*R*(1.0 - GG)
Q = R*pow(TT*(1.0 - alpha)/alpha, 1.0/(RR - 1.0))
a2eqss = a2max/(1.0 + Q*b)
weqss = fitness(a2eqss, a2eqss)
q2eqss = a2eqss*R2

for row, given, q2eqs, weqs in zip(axs, givens, q2eqss, weqss):
    budget0 = q2_budget*(1.0 - given)
    for ax, rho, ics, q2eq, weq in zip(row, rhos, icss, q2eqs, weqs):
        for ic in ics:
            ax.plot(q1_ic, ic, c='0.850')
        budget = budget0 + q2eq*given
        ax.plot(q1_budget, budget, c='green')
        ax.plot(q1_ic, icces(q1_ic, weq, rho), linewidth=2, c=cm.magma(weq))
        ax.set(xticks=[], yticks=[], xlim=(0.0, a1max*R1), ylim=(0.0, a2max*R2))
        ax.set_box_aspect(1)
axs[0, 0].set_title('a', fontsize=fslabel, weight='bold')
for ax, given in zip(axs[::every, 0], givens[::every]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)
for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
    ax.set_xlabel(round(log_es), fontsize=fstick)

# Indifference curves and discrete budget line

grid = outer_grid[0, 1].subgridspec(num, num, wspace=0, hspace=0)
axs = grid.subplots()

a2Ds = np.full([num, num], 0.0)
a2Cs = np.full([num, num], a2max/2.0)
Ts = fitness(a2Cs, a2Ds)
Rs = fitness(a2Cs, a2Cs)
Ps = fitness(a2Ds, a2Ds)
Ss = fitness(a2Ds, a2Cs)
xeqss = np.full([num, num], 0.0)
mask = (Ts < Rs) & (Ps < Ss)
xeqss[mask] = 1.0
mask = (Ts > Rs) & (Ps < Ss) & (Rs - Ss - Ts + Ps != 0.0)
xeqss[mask] = (Ps[mask] - Ss[mask])/(Rs[mask] - Ss[mask] - Ts[mask] + Ps[mask])
q2eqss = xeqss*R2*a2max/2.0
weqss = (Ts + Ss)*xeqss*(1.0 - xeqss) + Rs*xeqss*xeqss + Ps*(1.0 - xeqss)*(1.0 - xeqss)

for row, given, xeqs, q2eqs, weqs in zip(axs, givens, xeqss, q2eqss, weqss):
    budget0 = q2_budget*(1.0 - given)
    budget5 = budget0 + (a2max/2.0)*R2*given
    for ax, rho, ics, xeq, q2eq, weq in zip(row, rhos, icss, xeqs, q2eqs, weqs):
        for ic in ics:
            ax.plot(q1_ic, ic, c='0.850')
        ax.plot(q1_budget, budget0, c='green', linewidth=0.5, marker='o', alpha=1.0-xeq)
        ax.plot(q1_budget, budget5, c='green', linewidth=0.5, marker='o', alpha=xeq)
        ax.plot(q1_ic, icces(q1_ic, weq, rho), linewidth=2, c=cm.magma(weq))
        ax.set(xticks=[], yticks=[], xlim=(0.0, a1max*R1), ylim=(0.0, a2max*R2))
        ax.set_box_aspect(1)
axs[0, 0].set_title('b', fontsize=fslabel, weight='bold')
for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
    ax.set_xlabel(round(log_es), fontsize=fstick)

plt.savefig('icurves.png', dpi=100)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')