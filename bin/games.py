#! /usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.colors import ListedColormap
import numpy as np
import time

start_time = time.perf_counter ()

alpha = 0.50
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0
npoints = 128
if alpha == 0.50:
    vmax = 1.2
    limmatrix = a2max/2.0
else:
    vmax = 1.5
    limmatrix = a2max

num = 11    # Number of subplot rows and columns
every = int(num/2)
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

def fitness(x, y, given, rho):
    if isinstance(x, float): x = np.array([x])
    if isinstance(y, float): y = np.array([y])
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - given) + x*R2*given
    if rho == 0.0:
        w = q1*q2
        mask = (w > 0.0)
        w[mask] = pow(q1[mask], alpha)*pow(q2[mask], 1.0 - alpha)
    elif rho < 0.0:
        w = q1*q2
        mask = (w > 0.0)
        w[mask] = alpha*pow(q1[mask], rho) + (1.0 - alpha)*pow(q2[mask], rho)
        mask = (w > 0.0)
        w[mask] = pow(w[mask], 1.0/rho)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

b = a2max/a1max
givens = np.linspace(maxgiven, mingiven, num=num)
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)

# Figure properties

fslabel = 26 # Label font size
fstick = 18 # Tick font size

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(12, 12)) 
fig.supylabel("Partner's share of $\it{A}$", x=0.04, y=0.520, fontsize=fslabel)
fig.supxlabel("Substitutability of $\it{A}$", x=0.525, y=0.05, fontsize=fslabel)

outer_grid = fig.add_gridspec(2, 2, left=0.15, right=0.9, top=0.9, bottom=0.15)

# Continuous game types

grid = outer_grid[0, 0].subgridspec(num, num, wspace=0, hspace=0)
axs = grid.subplots()

a2 = np.linspace(0.0, limmatrix, num=npoints)
X, Y = np.meshgrid(a2, a2)
#grey = [0.4, 0.4, 0.4, 1.0]
#red = [1.0, 0.0, 0.0, 1.0]
#orange = [1.0, 165/265, 0.0, 1.0]
#cyan = [0.0, 1.0, 1.0, 1.0]
mycolors = ['cyan', 'blue', 'red', 'orange']
mycmap = ListedColormap(mycolors)

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        T = fitness(Y, X, given, rho)
        R = fitness(Y, Y, given, rho)
        P = fitness(X, X, given, rho)
        S = fitness(X, Y, given, rho)
        mask = (R < P)
        H = R[mask]
        R[mask] = P[mask]
        P[mask] = H
        H = T[mask]
        T[mask] = S[mask]
        S[mask] = H
        PD = ((T > R) & (P > S)).astype(int) 
        SD = ((T >= R) & (P <= S)).astype(int) 
        ND = ((T < R) & (P < S)).astype(int) 
        Z = PD*1 + SD*2 + ND*3 + 1
        Z = np.tril(Z, k=-1)
        Z = np.ma.masked_where(Z == 0.0, Z)
        Z = Z - 1
        cmap = cm.get_cmap(mycmap).copy()
        cmap.set_bad(color='white')
        ax.imshow(Z, origin='lower', cmap=cmap, vmin=0, vmax=3)
        #ax.set_facecolor('0.200')
        ax.set(xticks=[], yticks=[], xlim=(-9, npoints + 5), ylim=(-5, npoints + 9))

axs[0, 0].set_title('a', fontsize=fslabel, weight='bold')
for ax, given in zip(axs[::every, 0], givens[::every]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

# Continuous fitness curves

grid = outer_grid[1, 0].subgridspec(num, num, wspace=0, hspace=0)
axs = grid.subplots()

givens[0] = 0.999999
Rq = R2/R1
a2s = np.linspace(0.0, a2max, num=npoints)
RR, GG = np.meshgrid(rhos, givens)
TT = b*Rq*(1.0 - GG)
Q = Rq*pow(TT*(1.0 - alpha)/alpha, 1.0/(RR - 1.0))
a2eqss = a2max/(1.0 + Q*b)

for row, given, a2eqs in zip(axs, givens, a2eqss):
    for ax, rho, a2eq in zip(row, rhos, a2eqs):
        w = []
        for a2 in a2s:
            w.append(fitness(a2eq, a2, given, rho))
        weq = fitness(a2eq, a2eq, given, rho)
        ax.plot(a2s, w, linewidth=2, c=cm.magma(weq/vmax))
        ax.set(xticks=[], yticks=[], xlim=(0.0, a2max), ylim=(0.0, 2.0))
        #ax.set_facecolor('0.200')
        #ax.set_box_aspect(1)

axs[0, 0].set_title('c', fontsize=fslabel, weight='bold')
for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
    ax.set_xlabel(round(log_es), fontsize=fstick)
for ax, given in zip(axs[::every, 0], givens[::every]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

# Discrete

grid = outer_grid[0, 1].subgridspec(num, num, wspace=0, hspace=0)
axs0 = grid.subplots()
grid = outer_grid[1, 1].subgridspec(num, num, wspace=0, hspace=0)
axs1 = grid.subplots()

givens[0] = 1.0
a2 = np.linspace(0.0, a2max, num=3)
xaxis = [1, 2, 3, 4]

for row0, row1, given in zip(axs0, axs1, givens):
    for ax0, ax1, rho in zip(row0, row1, rhos):
        w = []
        T = fitness(a2[1], a2[0], given, rho)
        R = fitness(a2[1], a2[1], given, rho)
        P = fitness(a2[0], a2[0], given, rho)
        S = fitness(a2[0], a2[1], given, rho)
        if R < P:
            H = T
            T = S
            S = H
            H = R
            R = P
            P = H
        Su = fitness(a2[1], a2[2], given, rho)
        if Su > R:
            T = fitness(a2[2], a2[1], given, rho)
            R = fitness(a2[2], a2[2], given, rho)
            P = fitness(a2[1], a2[1], given, rho)
            S = Su
        if (T < R) & (P < S):
            x = 1.0
            weq = R
            rgb = 'orange'
        elif (T >= R) & (P <= S) & (R - S - T + P != 0.0):
            x = (P - S)/(R - S - T + P)
            weq = (T + S)*x*(1.0 - x) + R*x*x + P*(1.0 - x)*(1.0 - x)
            rgb = 'red'
        elif (T > R) & (P > S):
            x = 0.0
            weq = P
            rgb = 'blue'
        else:
            x = 0.0
            weq = P
            rgb = 'cyan'
        for a in a2:
            w.append(fitness(a2[0], a, given, rho)*(1.0 - x) + fitness(a2[1], a, given, rho)*x)   
        yaxis = [T, R, P, S]
        ax0.plot(xaxis, yaxis, color=rgb, marker='o', markerfacecolor='white', linewidth=2, markersize=3)
        ax0.set(xticks=[], yticks=[], xlim=(0, 5), ylim=(0.0, 2.0))
        #ax0.set_facecolor('0.200')
        #ax0.set_box_aspect(1)
        ax1.plot(a2, w, linewidth=2, c=cm.magma(weq/vmax))
        ax1.set(xticks=[], yticks=[], xlim=(0.0, a2max), ylim=(0.0, 2.0))
        #ax1.set_facecolor('0.200')
        #ax1.set_box_aspect(1)

axs0[0, 0].set_title('b', fontsize=fslabel, weight='bold')
axs1[0, 0].set_title('d', fontsize=fslabel, weight='bold')
for ax1, log_es in zip(axs1[-1, ::every], log_ess[::every]):
    ax1.set_xlabel(round(log_es), fontsize=fstick)

plt.savefig('games.png', dpi=200)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
