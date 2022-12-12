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

minalpha = 0.4
maxalpha = 0.6
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

num = 21    # Number of subplot rows and columns
filename = 'gamesdiscsalpha'
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

# Figure

traitlabels = ['Effort to get $\it{A}$', 'Help', 'Fitness']
traitvmaxs = [1.0, 2.0, 1.16]
fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

letters = [['a', 'b', 'c'],
            ['b', 'c', 'd'],
            ['e', 'f', 'g']]

def fitness(x, y):
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - given) + x*R2*given
    w = q1*q2
    mask = (w > 0.0) & (RR == 0.0)
    w[mask] = pow(q1[mask], AA[mask])*pow(q2[mask], 1.0 - AA[mask])
    mask = (w > 0.0) & (RR < 0.0)
    w[mask] = AA[mask]*pow(q1[mask], RR[mask]) + (1.0 - AA[mask])*pow(q2[mask], RR[mask])
    mask = (w > 0.0) & (RR < 0.0)
    w[mask] = pow(w[mask], 1.0/RR[mask])
    mask = (RR > 0.0)
    w[mask] = pow(AA[mask]*pow(q1[mask], RR[mask]) + (1.0 - AA[mask])*pow(q2[mask], RR[mask]), 1.0/RR[mask])
    return w

if mingiven != maxgiven:
    movie = True
    givens = np.linspace(maxgiven, mingiven, num=num)
    frames = []
else:
    movie = False 
    givens = np.array([mingiven])

alphas = np.linspace(maxalpha, minalpha, num=num)
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
b = a2max/a1max
RR, AA = np.meshgrid(rhos, alphas)
a20 = np.full([num, num], 0.0)
a21 = np.full([num, num], a2max/2.0)
a22 = np.full([num, num], a2max)
given = 0.0
w00 = fitness(a20, a20)
w11 = fitness(a21, a21)
w22 = fitness(a22, a22)
a2optimal = a20
mask = (w11 > w00)
a2optimal[mask] = a21[mask]
mask = (w22 > w11)
a2optimal[mask] = a22[mask]
woptimal = fitness(a2optimal, a2optimal)

minx = round(log_ess[0])
maxx = round(log_ess[-1])
miny = minalpha
maxy = maxalpha

xticklabels = [minx, round((minx + maxx)/2), maxx]
yticklabels = [miny, (miny + maxy)/2, maxy]
extent = 0, num, 0, num

for given in reversed(givens):

    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(18, 18))
    fig.delaxes(axs[0, 1])
    fig.delaxes(axs[0, 2])
    fig.supylabel("Partner's share of $\it{A}$", x=0.04, y=0.520, fontsize=fslabel)
    fig.supxlabel("Substitutability of $\it{A}$", x=0.525, y=0.05, fontsize=fslabel)

    if movie:
        fig.text(0.80, 0.80, f'given\n{round(given,2)}', fontsize=fstick+4, color='grey', ha='right')

    Z = np.full([num, num], 0.0)
    a20 = np.full([num, num], 0.0)
    a21 = np.full([num, num], a2max/2.0)
    a2 = w = x = a20
    w01 = fitness(a20, a21)
    w10 = fitness(a21, a20)
    w12 = fitness(a21, a22)
    w21 = fitness(a22, a21)
    w02 = fitness(a20, a22)
    w20 = fitness(a22, a20)

    mask = (woptimal == w00)
    T = w01
    R = w00
    P = w11
    S = w10
    x[mask] = (P[mask] - S[mask])/(R[mask] - S[mask] - T[mask] + P[mask])
    x[(x < 0.0)] = x[(x < 0.0)]*0.0
    x[(x > 1.0)] = pow(x[(x > 1.0)], 0.0)
    a2[mask] = a20[mask]*x[mask] + a21[mask]*(1.0 - x[mask])
    w[mask] = (T[mask] + S[mask])*x[mask]*(1.0 - x[mask]) + R[mask]*x[mask]*x[mask] + P[mask]*(1.0 - x[mask])*(1.0 - x[mask])
    Z[mask & (T < R) & (P < S)] = 0.9
    Z[mask & (T >= R) & (P <= S)] = 0.5
    Z[mask & (T > R) & (P > S)] = 0.1

    mask = (woptimal == w11) & (w20 > w22)
    T = w10
    R = w11
    P = w00
    S = w01
    x[mask] = (P[mask] - S[mask])/(R[mask] - S[mask] - T[mask] + P[mask])
    x[(x < 0.0)] = x[(x < 0.0)]*0.0
    x[(x > 1.0)] = pow(x[(x > 1.0)], 0.0)
    a2[mask] = a21[mask]*x[mask] + a20[mask]*(1.0 - x[mask])
    w[mask] = (T[mask] + S[mask])*x[mask]*(1.0 - x[mask]) + R[mask]*x[mask]*x[mask] + P[mask]*(1.0 - x[mask])*(1.0 - x[mask])
    Z[mask & (T < R) & (P < S)] = 0.9
    Z[mask & (T >= R) & (P <= S)] = 0.5
    Z[mask & (T > R) & (P > S)] = 0.1

    mask = (woptimal == w11) & (w20 <= w22)
    T = w12
    R = w11
    P = w22
    S = w21
    x[mask] = (P[mask] - S[mask])/(R[mask] - S[mask] - T[mask] + P[mask])
    x[(x < 0.0)] = x[(x < 0.0)]*0.0
    x[(x > 1.0)] = pow(x[(x > 1.0)], 0.0)
    a2[mask] = a21[mask]*x[mask] + a22[mask]*(1.0 - x[mask])
    w[mask] = (T[mask] + S[mask])*x[mask]*(1.0 - x[mask]) + R[mask]*x[mask]*x[mask] + P[mask]*(1.0 - x[mask])*(1.0 - x[mask])
    Z[mask & (T < R) & (P < S)] = 0.9
    Z[mask & (T >= R) & (P <= S)] = 0.5
    Z[mask & (T > R) & (P > S)] = 0.1

    mask = (woptimal == w22) & (w10 < w11)
    T = w21
    R = w22
    P = w11
    S = w12
    x[mask] = (P[mask] - S[mask])/(R[mask] - S[mask] - T[mask] + P[mask])
    x[(x < 0.0)] = x[(x < 0.0)]*0.0
    x[(x > 1.0)] = pow(x[(x > 1.0)], 0.0)
    a2[mask] = a22[mask]*x[mask] + a21[mask]*(1.0 - x[mask])
    w[mask] = (T[mask] + S[mask])*x[mask]*(1.0 - x[mask]) + R[mask]*x[mask]*x[mask] + P[mask]*(1.0 - x[mask])*(1.0 - x[mask])
    Z[mask & (T < R) & (P < S)] = 0.9
    Z[mask & (T >= R) & (P <= S)] = 0.5
    Z[mask & (T > R) & (P > S)] = 0.1

    mask = (woptimal == w22) & (w10 >= w11)
    T = w10
    R = w11
    P = w00
    S = w01
    x[mask] = (P[mask] - S[mask])/(R[mask] - S[mask] - T[mask] + P[mask])
    x[(x < 0.0)] = x[(x < 0.0)]*0.0
    x[(x > 1.0)] = pow(x[(x > 1.0)], 0.0)
    a2[mask] = a21[mask]*x[mask] + a20[mask]*(1.0 - x[mask])
    w[mask] = (T[mask] + S[mask])*x[mask]*(1.0 - x[mask]) + R[mask]*x[mask]*x[mask] + P[mask]*(1.0 - x[mask])*(1.0 - x[mask])
    Z[mask & (T < R) & (P < S)] = 0.9
    Z[mask & (T >= R) & (P <= S)] = 0.5
    Z[mask & (T > R) & (P > S)] = 0.1

    Mss = [[a2optimal, a2optimal*R2*given, woptimal], [a2, a2*R2*given, w]]

    for axrow, letterrow in zip(axs, letters):
        for ax, letter, traitlabel in zip(axrow, letterrow, traitlabels):
            ax.set(xticks=[0, num/2, num], yticks=[0, num/2, num], xticklabels=[], yticklabels=[])
            ax.text(0, num*1.035, letter, fontsize=fslabel, weight='bold')
            if ax.get_subplotspec().is_last_row():
                ax.set_xticklabels(xticklabels, fontsize=fstick)
            if ax.get_subplotspec().is_first_col():
                ax.set_yticklabels(yticklabels, fontsize=fstick) 

    axs[0, 0].imshow(Z, extent=extent, cmap='magma', vmin=0, vmax=1)

    for row, Ms in zip(axs[1:], Mss):
        for ax, M, traitvmax in zip(row, Ms, traitvmaxs):
            ax.imshow(M, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)

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
