#! /usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import time

start_time = time.perf_counter ()

alpha = 0.5
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0
npoints = 32
npoints_ic = 32

n_ic = 6    # Number of indifference curves

num = 11
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0
aC = 0.5
aD = 0.0

def a2maxw(x, given, rho):
    T = b*R*(1.0 - given)
    Q = R*pow(T*(1.0 - alpha)/alpha, 1.0/(rho - 1.0))
    a2 = (a2max - x*given*Q*b)/(1.0 + Q*b*(1.0 - given))
    a2 = 1.0 - a2
    return a2

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
givens[0] = 0.99
R = R2/R1
b = a2max/a1max

fslabel = 26 # Label font size
fstick = 18 # Tick font size

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(6, 18), constrained_layout=False) 
fig.supylabel("Partner's share of $\it{A}$", fontsize=fslabel)
fig.supxlabel("Substitutability of $\it{A}$", x=0.555, y=0.045, fontsize=fslabel)

outer_grid = fig.add_gridspec(3, 1, left=0.2, right=0.9)

# Plot indifference curves and budget line

top_grid = outer_grid[0].subgridspec(num, num, wspace=0, hspace=0)
axs = top_grid.subplots()

x_ics = np.linspace(0.0, R1*1.5, num=npoints_ic)

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        for w in np.linspace(0.5, 3.0, num=n_ic):
            ax.plot(x_ics, icces(w, rho), c='#dbdbdb')
        a2 = a2eq(given, rho)
        w = fitness(a2, a2, given, rho)
        ax.plot(x_ics, icces(w, rho), c='#7e7e7e')
        T = b*R*(1.0 - given)
        budget = R2*a2max*(1.0 - given) + a2*R2*given - T*x_ics
        ax.plot(x_ics, budget, c='orange')
        ax.set(xticks=[], yticks=[], xlim=(0, R1*1.5), ylim=(0, R2*1.5))
        ax.set_box_aspect(1)
axs[0, 0].set_title('a', fontsize=fslabel, weight='bold')
for ax, given in zip(axs[::5, 0], givens[::5]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

# Plot fitness landscapes

center_grid = outer_grid[1].subgridspec(num, num, wspace=0, hspace=0)
axs = center_grid.subplots()

x = np.linspace(0.001, 0.999, num=npoints)
a2partner = x
y = np.linspace(0.999, 0.001, num=npoints)
X, Y = np.meshgrid(x, y)

extent = 0, npoints, npoints, 0
for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        Z = fitness(X, Y, given, rho)
        Z_normed = Z/Z.max(axis=0)
        ax.imshow(Z_normed, extent=extent, cmap='magma', vmin=0, vmax=1.1)
        xaxis = a2partner*npoints
        yaxis = a2maxw(a2partner, given, rho)*npoints
        ax.plot(xaxis, yaxis, color='white')
        ax.set(xticks=[], yticks=[], xlim=(-0.5, npoints-0.5), ylim=(npoints-0.5, -0.5))
for ax, given in zip(axs[::5, 0], givens[::5]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)
axs[0, 0].set_title('b', fontsize=fslabel, weight='bold')

# Plot game types

bottom_grid = outer_grid[2].subgridspec(num, num, wspace=0, hspace=0)
axs = bottom_grid.subplots()

xaxis = [1, 2, 3, 4]
givens[0] = 1.0

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        R = wC1 = fitness(aC, aC, given, rho)
        P = wD0 = fitness(aD, aD, given, rho)
        T = wD1 = fitness(aC, aD, given, rho)
        S = wC0 = fitness(aD, aC, given, rho)
        yaxis = [T, R, P, S]
        if (T < R) and (P < S):
            rgb = 'orange'
        elif (T > R) and (P < S):
            rgb = 'red'
        elif (2.0*R > T + S):
            rgb = 'black'
        else:
            rgb = (1.0, 0.0, 1.0)
        ax.plot(xaxis, yaxis, color=rgb, marker='o', markerfacecolor='white', linewidth=1.0, markersize=3)
        ax.set(xticks=[], yticks=[], xlim=(0, 5), ylim=(0.0, 2.0))
        ax.set_box_aspect(1)

for ax, given in zip(axs[::5, 0], givens[::5]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)
for ax, log_es in zip(axs[-1, ::5], log_ess[::5]):
    ax.set_xlabel(round(log_es), fontsize=fstick)
axs[0, 0].set_title('c', fontsize=fslabel, weight='bold')

plt.savefig('games.png', dpi=100)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
