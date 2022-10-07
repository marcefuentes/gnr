#! /usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
import numpy as np
import time

start_time = time.perf_counter ()

alpha = 0.5
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0
npoints = 128
npoints_ic = 32

n_ic = 3    # Number of indifference curves

num = 11
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

def a2eq(given, rho):
    T = b*R*(1.0 - given)
    Q = R*pow(T*(1.0 - alpha)/alpha, 1.0/(rho - 1.0))
    a2 = a2max/(1.0 + Q*b)
    return a2

def fitness(x, y, given, rho):
    q1 = (1.0 - y)*R1
    q2 = y*R2*(1.0 - given) + x*R2*given
    if rho == 0.0:
        w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

def icces(w, rho):
    if rho == 0.0:
        ics = pow(w/pow(x_ics, alpha), 1.0/(1.0 - alpha))
    else:
        ics = pow((pow(w, rho) - alpha*pow(x_ics, rho))/(1.0 - alpha), 1.0/rho)
        if rho < 0.0:
            ics[pow(w, rho) <= alpha*pow(x_ics, rho)] = 1000.0
        else:
            ics[pow(w, rho) <= alpha*pow(x_ics, rho)] = -0.1
    return ics

log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
givens = np.linspace(maxgiven, mingiven, num=num)
givens[0] = 0.9999999
R = R2/R1
b = a2max/a1max
Ts = b*R*(1.0 - givens)
RR, TT = np.meshgrid(rhos, Ts)
Q = R*pow(TT*(1.0 - alpha)/alpha, 1.0/(RR - 1.0))

fslabel = 26 # Label font size
fstick = 18 # Tick font size

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(12, 12), constrained_layout=False) 
fig.supylabel("Partner's share of $\it{A}$", y=0.485, fontsize=fslabel)
fig.supxlabel("Substitutability of $\it{A}$", x=0.525, fontsize=fslabel)

outer_grid = fig.add_gridspec(2, 2, left=0.15, right=0.9, top=0.86, bottom=0.11)

# Top left: indifference curves and budget line

topleft_grid = outer_grid[0, 0].subgridspec(num, num, wspace=0, hspace=0)
axs = topleft_grid.subplots()

x_ics = np.linspace(0.0, R1*a1max, num=npoints_ic)
a2eqs = a2max/(1.0 + Q*b)
ws = np.linspace(0.5, 1.5, num=5)

for row, given, g in zip(axs, givens, a2eqs):
    for ax, rho, a2 in zip(row, rhos, g):
        for w in ws:
            ax.plot(x_ics, icces(w, rho), c='#dbdbdb')
        weq = fitness(a2, a2, given, rho)
        ax.plot(x_ics, icces(weq, rho), c=cm.magma(weq))
        T = b*R*(1.0 - given)
        budget = R2*a2max*(1.0 - given) + a2*R2*given - T*x_ics
        ax.plot(x_ics, budget, c='green', alpha=0.4)
        ax.set(xticks=[], yticks=[], xlim=(0, R1*a1max), ylim=(0, R2*a2max))
        ax.set_box_aspect(1)
axs[0, 0].set_title('a', fontsize=fslabel, weight='bold')
for ax, given in zip(axs[::5, 0], givens[::5]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

# Top right: fitness landscapes

topright_grid = outer_grid[0, 1].subgridspec(num, num, wspace=0, hspace=0)
axs = topright_grid.subplots()

a2 = np.linspace(0.001, 0.999, num=npoints)
X, Y = np.meshgrid(a2, a2)

for row, given, g in zip(axs, givens, Q):
    for ax, rho, q in zip(row, rhos, g):
        Z = fitness(X, Y, given, rho)
        Z_normed = Z/Z.max(axis=0)
        ax.imshow(Z_normed, cmap='magma', vmin=0, vmax=1.1)
        xaxis = npoints*a2-0.5
        a2maxw = npoints*(a2max - a2*given*q*b)/(1.0 + q*b*(1.0 - given))-0.5
        ax.plot(xaxis, a2maxw, color='white', alpha=0.7)
        ax.set(xticks=[], yticks=[], xlim=(-0.5, npoints-0.5), ylim=(-0.5, npoints-0.5))
axs[0, 0].set_title('b', fontsize=fslabel, weight='bold')

# Bottom left: extreme game types

bottomleft_grid = outer_grid[1, 0].subgridspec(num, num, wspace=0, hspace=0)
axs = bottomleft_grid.subplots()

xaxis = [1, 2, 3, 4]

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        #aC = a2eq(0.0, rho)
        #aD = a2eq(given, rho)
        aC = 0.5
        aD = 0.0
        R = fitness(aC, aC, given, rho)
        P = fitness(aD, aD, given, rho)
        T = fitness(aC, aD, given, rho)
        S = fitness(aD, aC, given, rho)
        yaxis = [T, R, P, S]
        if (T <= R) and (P <= S):
            rgb = 'orange'
        elif (T > R) and (P <= S):
            rgb = 'crimson'
        elif (2.0*R > T + S):
            rgb = (0.8-(R-P), 0.8-(R-P), 0.8-(R-P))
        else:
            rgb = 'cyan'
        ax.plot(xaxis, yaxis, color=rgb, marker='o', markerfacecolor='white', linewidth=1.0, markersize=3)
        ax.set(xticks=[], yticks=[], xlim=(0, 5), ylim=(0.0, 2.0))
        ax.set_box_aspect(1)

for ax, given in zip(axs[::5, 0], givens[::5]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)
for ax, log_es in zip(axs[-1, ::5], log_ess[::5]):
    ax.set_xlabel(round(log_es), fontsize=fstick)
axs[0, 0].set_title('c', fontsize=fslabel, weight='bold')

# Bottom right: many game types

bottomright_grid = outer_grid[1, 1].subgridspec(num, num, wspace=0, hspace=0)
axs = bottomright_grid.subplots()

a2 = np.linspace(0.000, 0.5, num=npoints)
X, Y = np.meshgrid(a2, a2)

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        R = fitness(Y, Y, given, rho)
        P = fitness(X, X, given, rho)
        T = fitness(Y, X, given, rho)
        S = fitness(X, Y, given, rho)
        TR = (T <= R).astype(int)
        PS = (P <= S).astype(int)
        Z = TR + PS + 1
        Z = np.tril(Z, k=-1)
        Z = np.ma.masked_where(Z == 0.0, Z)
        Z = Z - 1
        cmap = cm.get_cmap('inferno').copy()
        cmap.set_bad(color='white')
        ax.imshow(Z, origin='lower', cmap=cmap, vmin=0, vmax=2.2)
        ax.set(xticks=[], yticks=[], xlim=(-0.5, npoints-0.5), ylim=(-0.5, npoints-0.5))
for ax, log_es in zip(axs[-1, ::5], log_ess[::5]):
    ax.set_xlabel(round(log_es), fontsize=fstick)
axs[0, 0].set_title('d', fontsize=fslabel, weight='bold')

plt.savefig('games.png', dpi=100)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')