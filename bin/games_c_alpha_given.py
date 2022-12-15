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

minalpha = 0.1
maxalpha = 0.9
minlog_es = 1.0
maxlog_es = 1.0
mingiven = 0.0
maxgiven = 1.0

num = 21    # Number of subplot rows and columns
npoints = 128
filename = 'games_c_alpha_given'
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

def fitness(x, y, given, alpha, rho):
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - given) + x*R2*given
    w = q1*q2
    mask = (w > 0.0) & (rho == 0.0)
    w[mask] = pow(q1[mask], 1.0 - alpha[mask])*pow(q2[mask], alpha[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = (1.0 - alpha[mask])*pow(q1[mask], rho[mask]) + alpha[mask]*pow(q2[mask], rho[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = pow(w[mask], 1.0/rho[mask])
    mask = (rho > 0.0)
    w[mask] = pow((1.0 - alpha[mask])*pow(q1[mask], rho[mask]) + alpha[mask]*pow(q2[mask], rho[mask]), 1.0/rho[mask])
    return w

if minlog_es != maxlog_es:
    movie = True
    log_ess = np.linspace(minlog_es, maxlog_es, num=num)
    frames = []
else:
    movie = False 
    log_ess = np.array([minlog_es])

givens = np.linspace(maxgiven, mingiven, num=num)
givens[0] = 0.999999
alphas = np.linspace(minalpha, maxalpha, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
b = a2max/a1max
Rq = R2/R1
a2x = np.linspace(0.0, a2max, num=npoints)
a2y = np.linspace(a2max, 0.0, num=npoints)
X, Y = np.meshgrid(a2x, a2y)
Z = np.zeros([npoints, npoints])
AA, GG = np.meshgrid(alphas, givens)
AA, G0 = np.meshgrid(alphas, np.zeros([num]))
TT = b*Rq*(1.0 - GG)
T0 = b*Rq*(1.0 - G0)

minx = minalpha
maxx = maxalpha
miny = mingiven
maxy = maxgiven

xticklabels = [minx, (minx + maxx)/2, maxx]
yticklabels = [miny, (miny + maxy)/2, maxy]
extent = 0, num, 0, num
extentZ = 0, npoints*num, 0, npoints*num
black = np.full((npoints, npoints, 4), [0.2, 0.0, 0.2, 1.0])
cyan = np.full((npoints, npoints, 4), [0.0, 1.0, 1.0, 1.0])
white = np.full((npoints, npoints, 4), [1.0, 1.0, 1.0, 1.0])
green = np.full((npoints, npoints, 4), [0.0, 1.0, 0.0, 1.0])

for rho in rhos:

    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(18, 18))
    fig.delaxes(axs[0, 1])
    fig.delaxes(axs[0, 2])
    fig.supylabel("Partner's share of $\it{A}$", x=0.04, y=0.520, fontsize=fslabel)
    fig.supxlabel("Value of $\it{A}$", x=0.525, y=0.05, fontsize=fslabel)

    if movie:
        fig.text(0.80, 0.80, f'log(es)\n{-log(1.0-rho, 2):4.1}', fontsize=fstick+4, color='grey', ha='right')

    Q0 = Rq*pow(T0*AA/(1.0 - AA), 1.0/(rho - 1.0))
    a20ss = a2max/(1.0 + Q0*b)
    Q = Rq*pow(TT*AA/(1.0 - AA), 1.0/(rho - 1.0))
    a2eqss = a2max/(1.0 + Q*b)

    Rh = np.full([npoints, npoints], rho)
    Zss = np.empty((0, npoints*num, 4))
    for given in givens:
        G = np.full([npoints, npoints], given)
        Zs = np.empty((npoints, 0, 4))
        for alpha in alphas:
            Z = np.full([npoints, npoints, 4], [0.0, 1.0, 0.0, 1.0])
            A = np.full([npoints, npoints], alpha)
            T = fitness(Y, X, G, A, Rh)
            R = fitness(Y, Y, G, A, Rh)
            P = fitness(X, X, G, A, Rh)
            S = fitness(X, Y, G, A, Rh)
            mask = (R < P)
            H = R[mask]
            R[mask] = P[mask]
            P[mask] = H
            H = T[mask]
            T[mask] = S[mask]
            S[mask] = H
            mask = (T > R) & (P > S)
            Z[mask] = black[mask]
            mask = (T >= R) & (P <= S)
            Z[mask] = cyan[mask]
            mask = (T < R) & (P < S)
            Z[mask] = white[mask]
            #Z = np.tril(Z, k=-1)
            Z = np.ma.masked_where(Z == 0.0, Z)
            Zs = np.append(Zs, Z, axis=1)
        Zss = np.append(Zss, Zs, axis=0)

    Rh = np.full([num, num], rho)
    Mss = [[a20ss, a20ss*R2*GG, fitness(a20ss, a20ss, GG, AA, Rh)], [a2eqss, a2eqss*R2*GG, fitness(a2eqss, a2eqss, GG, AA, Rh)]]

    for axrow, letterrow in zip(axs, letters):
        for ax, letter, traitlabel in zip(axrow, letterrow, traitlabels):
            if ax.get_subplotspec().is_first_row():
                ax.set(xticks=[0, npoints*num/2, npoints*num], yticks=[0, npoints*num/2, npoints*num], xticklabels=[], yticklabels=yticklabels)
                #ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)
                ax.text(0, npoints*num*1.035, letter, fontsize=fslabel, weight='bold')
            else:
                ax.set(xticks=[0, num/2, num], yticks=[0, num/2, num], xticklabels=[], yticklabels=[])
                ax.text(0, num*1.035, letter, fontsize=fslabel, weight='bold')
            if ax.get_subplotspec().is_last_row():
                ax.set_xticklabels(xticklabels, fontsize=fstick)
            if ax.get_subplotspec().is_first_col():
                ax.set_yticklabels(yticklabels, fontsize=fstick) 

    axs[0, 0].imshow(Zss, extent=extentZ)

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
