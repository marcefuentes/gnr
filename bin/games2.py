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
vmax = 1.5

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
        w[w>0.0] = pow(q1[w>0.0], alpha)*pow(q2[w>0.0], 1.0 - alpha)
    elif rho < 0.0:
        w = q1*q2
        w[w>0.0] = alpha*pow(q1[w>0.0], rho) + (1.0 - alpha)*pow(q2[w>0.0], rho)
        w[w>0.0] = pow(w[w>0.0], 1.0/rho)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

b = a2max/a1max
givens = np.linspace(maxgiven, mingiven, num=num)
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)

fslabel = 26 # Label font size
fstick = 18 # Tick font size

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(12, 12), constrained_layout=False) 
fig.supylabel("Partner's share of $\it{A}$", x=0.04, y=0.520, fontsize=fslabel)
fig.supxlabel("Substitutability of $\it{A}$", x=0.525, y=0.05, fontsize=fslabel)

outer_grid = fig.add_gridspec(2, 2, left=0.15, right=0.9, top=0.9, bottom=0.15)

# Continuous game types

grid = outer_grid[0, 0].subgridspec(num, num, wspace=0, hspace=0)
axs = grid.subplots()

a2 = np.linspace(0.000, a2max, num=npoints)
X, Y = np.meshgrid(a2, a2)
grey = [0.4, 0.4, 0.4, 1.0]
red = [1.0, 0.0, 0.0, 1.0]
orange = [1.0, 165/265, 0.0, 1.0]
cyan = [0.0, 1.0, 1.0, 1.0]
mycolors = [cyan, grey, red, orange]
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
        SD = ((T > R) & (P < S)).astype(int) 
        ND = ((T < R) & (P < S)).astype(int) 
        Z = PD*1 + SD*2 + ND*3 + 1
        Z = np.tril(Z, k=-1)
        Z = np.ma.masked_where(Z == 0.0, Z)
        Z = Z - 1
        cmap = cm.get_cmap(mycmap).copy()
        cmap.set_bad(color='white')
        ax.imshow(Z, origin='lower', cmap=cmap, vmin=0, vmax=3)
        ax.set(xticks=[], yticks=[], xlim=(-9, npoints + 5), ylim=(-5, npoints + 9))
axs[0, 0].set_title('a', fontsize=fslabel, weight='bold')
for ax, given in zip(axs[::every, 0], givens[::every]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

# Discrete game types

grid = outer_grid[0, 1].subgridspec(num, num, wspace=0, hspace=0)
axs = grid.subplots()

a2 = np.array([0.0, a2max/2.0, a2max])
X, Y = np.meshgrid(a2, a2)
print(X, Y)
xaxis = [1, 2, 3, 4]

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        Z = fitness(X, Y, given, rho)
        R = max(Z[0, 0], Z[1, 1])
        if R == Z[1, 1]:
            T = Z[0, 1]
            S = Z[1, 0]
            P = Z[0, 0]
        else:
            T = Z[1, 0]
            S = Z[0, 1]
            P = Z[1, 1]
        yaxis = [T, R, P, S]
        if (T < R) and (P < S):
            rgb = 'orange'
        elif (T > R) and (P < S):
            rgb = 'red'
        elif (T > R) and (P > S) and (2.0*R > T + S):
            rgb = 'black'
        else:
            rgb = 'cyan'
        ax.plot(xaxis, yaxis, color=rgb, marker='o', markerfacecolor='white', linewidth=2, markersize=3)
        R = max(Z[1, 1], Z[2, 2])
        if R == Z[2, 2]:
            T = Z[1, 2]
            S = Z[2, 1]
            P = Z[1, 1]
        else:
            T = Z[2, 1]
            S = Z[1, 2]
            P = Z[2, 2]
        yaxis = [T, R, P, S]
        if (T < R) and (P < S):
            rgb = 'orange'
        elif (T > R) and (P < S):
            rgb = 'red'
        elif (T > R) and (P > S) and (2.0*R > T + S):
            rgb = 'black'
        else:
            rgb = 'cyan'
        ax.plot(xaxis, yaxis, color=rgb, marker='o', markerfacecolor='white', linewidth=2, markersize=3)
        ax.set(xticks=[], yticks=[], xlim=(0, 5), ylim=(0.0, 2.0))
        ax.set_box_aspect(1)
axs[0, 0].set_title('b', fontsize=fslabel, weight='bold')

# Continuous

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
        ax.set_facecolor('0.200')
        ax.set_box_aspect(1)
axs[0, 0].set_title('c', fontsize=fslabel, weight='bold')
for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
    ax.set_xlabel(round(log_es), fontsize=fstick)
for ax, given in zip(axs[::every, 0], givens[::every]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

givens[0] = 1.0

# Discrete

grid = outer_grid[1, 1].subgridspec(num, num, wspace=0, hspace=0)
axs = grid.subplots()

a2 = np.linspace(0.0, a2max, num=3)

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        w = []
        T = fitness(a2[1], a2[0], given, rho)
        R = fitness(a2[1], a2[1], given, rho)
        P = fitness(a2[0], a2[0], given, rho)
        S = fitness(a2[0], a2[1], given, rho)
        if (T < R) & (P < S):
            Tu = fitness(a2[2], a2[1], given, rho)
            Ru = fitness(a2[2], a2[2], given, rho)
            Pu = fitness(a2[1], a2[1], given, rho)
            Su = fitness(a2[1], a2[2], given, rho)
            if (Tu < Ru) & (Pu < Su):
                x = 1.0
                weq = Ru
            elif (Tu > Ru) & (Pu < Su) & (Ru - Su - Tu + Pu != 0.0):
                x = (Pu - Su)/(Ru - Su - Tu + Pu)
                weq = (Tu + Su)*x*(1.0 - x) + Ru*x*x + Pu*(1.0 - x)*(1.0 - x)
            else:
                x = 0.0
                weq = Pu
            for a in a2:
                w.append(fitness(a2[1], a, given, rho)*(1.0 - x) + fitness(a2[2], a, given, rho)*x)   
        else:
            if (T > R) & (P < S) & (R - S - T + P != 0.0):
                x = (P - S)/(R - S - T + P)
                weq = (T + S)*x*(1.0 - x) + R*x*x + P*(1.0 - x)*(1.0 - x)
            else:
                x = 0.0
                weq = P
            for a in a2:
                w.append(fitness(a2[0], a, given, rho)*(1.0 - x) + fitness(a2[1], a, given, rho)*x)   
        ax.plot(a2, w, linewidth=2, c=cm.magma(weq/vmax))
        ax.set(xticks=[], yticks=[], xlim=(0.0, a2max), ylim=(0.0, 2.0))
        ax.set_facecolor('0.200')
        ax.set_box_aspect(1)
axs[0, 0].set_title('d', fontsize=fslabel, weight='bold')
for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
    ax.set_xlabel(round(log_es), fontsize=fstick)

plt.savefig('games.png', dpi=200)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
