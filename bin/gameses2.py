#! /usr/bin/env python

from glob import glob
import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.colors import ListedColormap
import numpy as np
import time

start_time = time.perf_counter ()

filename = 'gameses2'
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0
npoints = 128
#vmax = 1.2
vmax = 1.5
num = 11    # Number of subplot rows and columns
every = int(num/2)
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

# Figure

fslabel = 26 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

letters = [['a', 'b'],
            ['c', 'd'],
            ['e', 'f'],
            ['g', 'h']]

movie = False
if movie:
    alphas = np.linspace(0.1, 0.9, num=11)
    frames = []
else:
    alphas = np.linspace(0.5, 0.5, num=1)

def fitness(x, y):
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

# Continuous

givens[0] = 0.999999
Rq = R2/R1
a2x = np.linspace(0.0, a2max, num=npoints)
a2y = np.linspace(a2max, 0.0, num=npoints)
X, Y = np.meshgrid(a2x, a2y)
Z = np.full([npoints, npoints], 0.0)
RR, GG = np.meshgrid(rhos, givens)
RR, G0 = np.meshgrid(rhos, np.full(num, 0.0))
TT = b*Rq*(1.0 - GG)
T0 = b*Rq*(1.0 - G0)

extent0 = 0, num, 0, num
extent = 0, npoints, 0, npoints

# Discrete

a2 = np.linspace(0.0, a2max, num=3)
xaxis = [1, 2, 3, 4]

for alpha in alphas:

    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 18))
    fig.supylabel("Partner's share of $\it{A}$", x=0.04, y=0.520, fontsize=fslabel)
    fig.supxlabel("Substitutability of $\it{A}$", x=0.525, y=0.05, fontsize=fslabel)

    if movie: fig.text(0.93, 0.02, f'alpha = {alpha}', fontsize=fstick, color='grey', ha='right')

    # Continuous

    Q = Rq*pow(TT*(1.0 - alpha)/alpha, 1.0/(RR - 1.0))
    a2eqss = a2max/(1.0 + Q*b)
    Q0 = Rq*pow(T0*(1.0 - alpha)/alpha, 1.0/(RR - 1.0))
    a20ss = a2max/(1.0 + Q0*b)

    Zss = np.empty((0, npoints*num))
    for given, a2eqs in zip(givens, a2eqss):
        Zs = np.empty((npoints, 0))
        for rho, a2eq in zip(rhos, a2eqs):
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
            Z[(T > R) & (P > S)] = 0.1
            Z[(T >= R) & (P <= S)] = 0.5
            Z[(T < R) & (P < S)] = 0.9
            Z[(T == R) & (P == R)] = 1.0
            #Z = np.tril(Z, k=-1)
            Z = np.ma.masked_where(Z == 0.0, Z)
            Zs = np.append(Zs, Z, axis=1)
        Zss = np.append(Zss, Zs, axis=0)

    axs[0, 0].imshow(Zss, origin='lower', cmap='magma', vmin=0, vmax=1)
    axs[0, 0].set(xticks=[], yticks=[], xlim=(0, npoints*num), ylim=(npoints*num, 0))
    axs[1, 0].imshow(a2eqss, origin='lower', extent=extent0, cmap='magma', vmin=0, vmax=1)
    axs[1, 0].set(xticks=[], yticks=[], xlim=(0, num), ylim=(num, 0))
    axs[2, 0].imshow(a20ss, origin='lower', extent=extent0, cmap='magma', vmin=0, vmax=1)
    axs[2, 0].set(xticks=[], yticks=[], xlim=(0, num), ylim=(num, 0))

    # Discrete

    givens[0] = 1.0


    if movie:
        plt.savefig('temp.png', transparent=False)
        frames.append(iio.imread('temp.png'))
        os.remove('temp.png')
    else:
        plt.savefig(filename + '.png', transparent=False)
if movie:
    iio.mimsave(filename + '.gif', frames)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
