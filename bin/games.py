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

log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
givens = np.linspace(maxgiven, mingiven, num=num)
R = R2/R1
b = a2max/a1max
givens[0] = 0.9999999
Ts = b*R*(1.0 - givens)
RR, TT = np.meshgrid(rhos, Ts)
Q = R*pow(TT*(1.0 - alpha)/alpha, 1.0/(RR - 1.0))


fslabel = 26 # Label font size
fstick = 18 # Tick font size

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(12, 6), constrained_layout=False) 
fig.supylabel("Partner's share of $\it{A}$", y=0.520, fontsize=fslabel)
fig.supxlabel("Substitutability of $\it{A}$", x=0.525, fontsize=fslabel)

outer_grid = fig.add_gridspec(1, 2, left=0.15, right=0.9, top=0.86, bottom=0.176)

# Fitness landscapes

grid = outer_grid[0, 0].subgridspec(num, num, wspace=0, hspace=0)
axs = grid.subplots()

a2 = np.linspace(0.0, a2max, num=npoints)
X, Y = np.meshgrid(a2, a2)
givens[0] = 0.9999999

for row, given, g in zip(axs, givens, Q):
    for ax, rho, q in zip(row, rhos, g):
        Z = fitness(X, Y, given, rho)
        Z_normed = Z/Z.max(axis=0)
        ax.imshow(Z_normed, origin='lower', cmap='Greys_r', vmin=0, vmax=1.0)
        a2maxw = (a2max - a2*given*q*b)/(1.0 + q*b*(1.0 - given))
        ax.plot(a2*npoints, a2maxw*npoints, color='orange')
        ax.set(xticks=[], yticks=[], xlim=(-0.5, npoints-0.5), ylim=(-0.5, npoints-0.5))
for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
    ax.set_xlabel(round(log_es), fontsize=fstick)
axs[0, 0].set_title('a', fontsize=fslabel, weight='bold')
for ax, given in zip(axs[::every, 0], givens[::every]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

# Discrete game types

grid = outer_grid[0, 1].subgridspec(num, num, wspace=0, hspace=0)
axs = grid.subplots()

a2 = np.array([a2max/2.0, 0.0])
X, Y = np.meshgrid(a2, a2)
givens[0] = 1.0
xaxis = [1, 2, 3, 4]

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        Z = fitness(X, Y, given, rho)
        T = Z[1, 0]
        R = Z[0, 0]
        P = Z[1, 1]
        S = Z[0, 1]
        yaxis = [T, R, P, S]
        if (T < R) and (P < S):
            rgb = 'orange'
        elif (T > R) and (P < S):
            rgb = 'red'
        elif (T > R) and (P > S) and (2.0*R > T + S):
            rgb = 'black'
        else:
            rgb = 'cyan'
        ax.plot(xaxis, yaxis, color=rgb, marker='o', markerfacecolor='white', linewidth=1.0, markersize=3)
        ax.set(xticks=[], yticks=[], xlim=(0, 5), ylim=(0.0, 2.0))
        ax.set_box_aspect(1)
for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
    ax.set_xlabel(round(log_es), fontsize=fstick)
axs[0, 0].set_title('b', fontsize=fslabel, weight='bold')

# Bottom right: many game types

#bottomright_grid = outer_grid[1, 1].subgridspec(num, num, wspace=0, hspace=0)
#axs = bottomright_grid.subplots()

#a2 = np.linspace(0.000, 0.5, num=npoints)
#X, Y = np.meshgrid(a2, a2)
#mycolors = [[0.4, 0.4, 0.4, 1.0], [1.0, 0.0, 0.0, 1.0], [1.0, 165/265, 0.0, 1.0]]
#mycmap = ListedColormap(mycolors)

#for row, given in zip(axs, givens):
#    for ax, rho in zip(row, rhos):
#        R = fitness(Y, Y, given, rho)
#        P = fitness(X, X, given, rho)
#        T = fitness(Y, X, given, rho)
#        S = fitness(X, Y, given, rho)
#        TR = (T <= R).astype(int)
#        PS = (P <= S).astype(int)
#        Z = TR + PS + 1
#        Z = np.tril(Z, k=-1)
#        Z = np.ma.masked_where(Z == 0.0, Z)
#        Z = Z - 1
#        cmap = cm.get_cmap(mycmap).copy()
#        cmap.set_bad(color='white')
#        ax.imshow(Z, origin='lower', cmap=cmap, vmin=0, vmax=2)
#        ax.set(xticks=[], yticks=[], xlim=(-0.5, npoints/2.0 - 0.5), ylim=(-0.5, npoints/2.0 - 0.5))
#for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
#    ax.set_xlabel(round(log_es), fontsize=fstick)
#axs[0, 0].set_title('d', fontsize=fslabel, weight='bold')

plt.savefig('games.png', dpi=100)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
