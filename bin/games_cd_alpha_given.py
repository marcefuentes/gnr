#! /usr/bin/env python

from glob import glob
from math import log
from matplotlib import cm
from matplotlib.colors import ListedColormap
import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import time

start_time = time.perf_counter ()

minalpha = 0.4
maxalpha = 0.6
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

num = 11    # Number of subplot rows and columns
every = int(num/2)
npoints = 128
filename = 'games_cd_alpha_given'
vmax = 1.5
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

# Figure

fslabel = 26 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def fitness(x, y):
    if isinstance(x, float): x = np.array([x])
    if isinstance(y, float): y = np.array([y])
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - given) + x*R2*given
    if rho == 0.0:
        w = q1*q2
        mask = (w > 0.0)
        w[mask] = pow(q1[mask], 1.0 - alpha)*pow(q2[mask], alpha)
    elif rho < 0.0:
        w = q1*q2
        mask = (w > 0.0)
        w[mask] = (1.0 - alpha)*pow(q1[mask], rho) + alpha*pow(q2[mask], rho)
        mask = (w > 0.0)
        w[mask] = pow(w[mask], 1.0/rho)
    else:
        w = pow((1.0 - alpha)*pow(q1, rho) + alpha*pow(q2, rho), 1.0/rho)
    return w

if minlog_es != maxlog_es:
    movie = True
    log_ess = np.linspace(minlog_es, maxlog_es, num=11)
    frames = []
else:
    movie = False
    log_ess = np.linspace(minlog_es, minlog_es, num=1)

alphas = np.linspace(maxalpha, minalpha, num=num)
givens = np.linspace(maxgiven, mingiven, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
b = a2max/a1max

# Continuous

givens[0] = 0.999999
Rq = R2/R1
a2s = np.linspace(0.0, a2max, num=npoints)
X, Y = np.meshgrid(a2s, a2s)
Z = np.full([npoints, npoints], 0.0)
AA, GG = np.meshgrid(alphas, givens)
TT = b*Rq*(1.0 - GG)
Q = np.full([num, num], 0.0)
a2eqss = np.full([num, num], 0.0)

extent = 0, npoints, 0, npoints

# Discrete

a2 = np.linspace(0.0, a2max, num=3)
xaxis = [1, 2, 3, 4]

for rho in rhos:

    fig = plt.figure(figsize=(12, 12)) 
    fig.supylabel("Partner's share of $\it{A}$", x=0.04, y=0.520, fontsize=fslabel)
    fig.supxlabel("Value of $\it{A}$", x=0.525, y=0.05, fontsize=fslabel)

    outer_grid = fig.add_gridspec(2, 2, left=0.15, right=0.9, top=0.9, bottom=0.15)

    if movie: fig.text(0.93, 0.02, f'log_es = {-log(1-rho)/log(2)}', fontsize=fstick, color='grey', ha='right')

    # Continuous

    grid = outer_grid[0, 0].subgridspec(num, num, wspace=0, hspace=0)
    ax0s = grid.subplots()
    grid = outer_grid[1, 0].subgridspec(num, num, wspace=0, hspace=0)
    ax1s = grid.subplots()

    Q = Rq*pow(TT*alpha/(1.0 - alpha), 1.0/(RR - 1.0))
    a2eqss = a2max/(1.0 + Q*b)

    for row0, row1, given, a2eqs in zip(ax0s, ax1s, givens, a2eqss):
        for ax0, ax1, alpha, a2eq in zip(row0, row1, alphas, a2eqs):
            T = fitness(Y, X)
            R = fitness(Y, Y)
            P = fitness(X, X)
            S = fitness(X, Y)
            mask = (R < P)
            H = R[mask]
            R[mask] = P[mask]
            P[mask] = H
            H = T[mask]
            T[mask] = S[mask]
            S[mask] = H
            Z[(T > R) and (P > S)] = 0.1
            Z[(T >= R) and (P <= S)] = 0.5
            Z[(T < R) and (P < S)] = 0.9
            Z[(T == R) and (P == R)] = 1.0
            #Z = np.tril(Z, k=-1)
            Z = np.ma.masked_where(Z == 0.0, Z)
            ax0.imshow(Z, origin='lower', extent=extent, cmap='magma', vmin=0, vmax=1)
            ax0.set(xticks=[], yticks=[], xlim=(-11, npoints + 11), ylim=(-11, npoints + 11))
            w = fitness(a2eq, a2s)
            weq = fitness(a2eq, a2eq)
            ax1.plot(a2s, w, linewidth=3, c=cm.magma(weq/vmax))
            ax1.set(xticks=[], yticks=[], xlim=(0.0, a2max), ylim=(0.0, 2.0))

    ax0s[0, 0].set_title('a', fontsize=fslabel, weight='bold')
    for ax0, given in zip(ax0s[::every, 0], givens[::every]):
        ax0.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)
    ax1s[0, 0].set_title('c', fontsize=fslabel, weight='bold')
    for ax1, alpha in zip(ax1s[-1, ::every], alphas[::every]):
        ax1.set_xlabel(round(alpha,1), fontsize=fstick)
    for ax1, given in zip(ax1s[::every, 0], givens[::every]):
        ax1.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

    # Discrete

    grid = outer_grid[0, 1].subgridspec(num, num, wspace=0, hspace=0)
    ax0s = grid.subplots()
    grid = outer_grid[1, 1].subgridspec(num, num, wspace=0, hspace=0)
    ax1s = grid.subplots()

    givens[0] = 1.0

    for row0, row1, given in zip(ax0s, ax1s, givens):
        for ax0, ax1, alpha in zip(row0, row1, alphas):
            a2low = a2[0]
            T = fitness(a2[1], a2[0])
            R = fitness(a2[1], a2[1])
            P = fitness(a2[0], a2[0])
            S = fitness(a2[0], a2[1])
            if R < P:
                H = T
                T = S
                S = H
                H = R
                R = P
                P = H
            Su = fitness(a2[1], a2[2])
            if Su > R:
                a2low = a2[1]
                T = fitness(a2[2], a2[1])
                R = fitness(a2[2], a2[2])
                P = fitness(a2[1], a2[1])
                S = Su
            if (T < R) and (P < S):
                x = 1.0
                weq = R
                rgb = 0.9
            elif (T >= R) and (P <= S) and (R - S - T + P != 0.0):
                x = (P - S)/(R - S - T + P)
                weq = (T + S)*x*(1.0 - x) + R*x*x + P*(1.0 - x)*(1.0 - x)
                rgb = 0.5
            elif (T > R) and (P > S):
                x = 0.0
                weq = P
                rgb = 0.1
            else:
                x = 0.0
                weq = P
                rgb = 0.0
            w = fitness(a2low, a2)*(1.0 - x) + fitness(a2low + a2max/2.0, a2)*x   
            yaxis = [T, R, P, S]
            ax0.plot(xaxis, yaxis, c=cm.magma(rgb), marker='o', markerfacecolor='white', linewidth=3, markersize=3)
            ax0.set(xticks=[], yticks=[], xlim=(0, 5), ylim=(0.0, 2.0))
            ax1.plot(a2, w, c=cm.magma(weq/vmax), linewidth=3)
            ax1.set(xticks=[], yticks=[], xlim=(0.0, a2max), ylim=(0.0, 2.0))

    ax0s[0, 0].set_title('b', fontsize=fslabel, weight='bold')
    ax1s[0, 0].set_title('d', fontsize=fslabel, weight='bold')
    for ax1, alpha in zip(ax1s[-1, ::every], alphas[::every]):
        ax1.set_xlabel(round(alpha, 1), fontsize=fstick)

    if movie:
        plt.savefig('temp.png', transparent=False)
        frames.append(iio.imread('temp.png'))
        os.remove('temp.png')
    else:
        plt.savefig(filename + '.png', transparent=False)

    plt.close()

if movie:
    iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
