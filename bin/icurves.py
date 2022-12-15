#! /usr/bin/env python

from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import time

start_time = time.perf_counter ()

npoints = 128
alpha = 0.5
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

num = 5    # Number of subplot rows and columns
#markersize = 10
n_ic = 5    # Number of indifference curves
every = int(num/2)
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

def fitness(x, y):
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - GG) + x*R2*GG
    w = q1*q2
    mask = (w > 0.0) & (RR == 0.0)
    w[mask] = pow(q1[mask], 1.0 - alpha)*pow(q2[mask], alpha)
    mask = (w > 0.0) & (RR < 0.0)
    w[mask] = (1.0 - alpha)*pow(q1[mask], RR[mask]) + alpha*pow(q2[mask], RR[mask])
    mask = (w > 0.0) & (RR < 0.0)
    w[mask] = pow(w[mask], 1.0/RR[mask])
    mask = (RR > 0.0)
    w[mask] = pow((1.0 - alpha)*pow(q1[mask], RR[mask]) + alpha*pow(q2[mask], RR[mask]), 1.0/RR[mask])
    return w

def icces(q, w):
    if rho == 0.0:
        q2 = np.piecewise(q, [q == 0.0, q > 0.0], [1000.0, lambda i: pow(w/pow(i, 1.0 - alpha), 1.0/alpha)])
    elif rho < 0.0:
        q2 = np.piecewise(q, [q == 0.0, q > 0.0], [1000.0, lambda i: np.piecewise(i, [pow(w, rho) <= (1.0 - alpha)*pow(i, rho), pow(w, rho) > (1.0 - alpha)*pow(i, rho)], [1000.0, lambda j: pow((pow(w, rho) - (1.0 - alpha)*pow(j, rho))/alpha, 1.0/rho)])])
    else:
        q2 = np.piecewise(q, [pow(w, rho) <= (1.0 - alpha)*pow(q, rho), pow(w, rho) > (1.0 - alpha)*pow(q, rho)], [-0.1, lambda i: pow((pow(w, rho) - (1.0 - alpha)*pow(i, rho))/alpha, 1.0/rho)])
    return q2

b = a2max/a1max
a1_budget = np.linspace(0.0, a1max, num=3)
q2_budget = (a2max - b*a1_budget)*R2
q1_budget = a1_budget*R1
q1_ic = np.linspace(0.0, a1max*R1, num=npoints)
givens = np.linspace(maxgiven, mingiven, num=num)
givens[0] = 0.9999999
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
RR, GG = np.meshgrid(rhos, givens)
wis = np.linspace(2.0/(n_ic + 1), 2.0*n_ic/(n_ic + 1), num=n_ic)
icss = []
for rho in rhos:
    ics = []
    for w in wis:
        ics.append(icces(q1_ic, w))
    icss.append(ics)

fslabel = 26 # Label font size
fstick = 18 # Tick font size

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(6, 6)) 
fig.supylabel("Partner's share of $\it{A}$", x=0.05, y=0.52, fontsize=fslabel)
fig.supxlabel("Substitutability of $\it{A}$", x=0.555, fontsize=fslabel)

grid = fig.add_gridspec(nrows=num, ncols=num, left=0.22, right=0.9, top=0.86, bottom=0.176, wspace=0, hspace=0)

axs = grid.subplots()

Rq = R2/R1
TT = b*Rq*(1.0 - GG)
Q = Rq*pow(TT*alpha/(1.0 - alpha), 1.0/(RR - 1.0))
a2ss = a2max/(1.0 + Q*b)
wss = fitness(a2ss, a2ss)
q2ss = a2ss*R2

for row, given, q2s, ws in zip(axs, givens, q2ss, wss):
    budget0 = q2_budget*(1.0 - given)
    for ax, rho, ics, q2eq, weq in zip(row, rhos, icss, q2s, ws):
        for ic in ics:
            ax.plot(q1_ic, ic, c='0.850')
        budget = budget0 + q2eq*given
        ax.plot(q1_budget, budget, c='green')
        ax.plot(q1_ic, icces(q1_ic, weq), linewidth=2, c=cm.magma(weq))
        ax.set(xticks=[], yticks=[], xlim=(0.0, a1max*R1), ylim=(0.0, a2max*R2))
        #ax.set_box_aspect(1)
for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
    ax.set_xlabel(round(log_es), fontsize=fstick)
for ax, given in zip(axs[::every, 0], givens[::every]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

plt.savefig('icurves.png', dpi=100)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
