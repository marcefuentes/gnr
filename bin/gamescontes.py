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

filename = 'gamescontes'
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

traitlabels = ['Effort to get $\it{A}$', 'Help', 'Fitness']
traitvmaxs = [1.0, 2.0, 2.0]
fslabel = 26 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

letters = [['a', 'b', 'c'],
            ['d', 'e', 'f'],
            ['g', 'h', 'i']]

movie = False
if movie:
    alphas = np.linspace(0.1, 0.9, num=11)
    frames = []
else:
    alphas = np.linspace(0.5, 0.5, num=1)

def fitness(x, y, given, alpha, rho):
    if isinstance(x, float): x = np.array([x])
    if isinstance(y, float): y = np.array([y])
    if isinstance(given, float): given = np.array([given])
    if isinstance(alpha, float): alpha = np.array([alpha])
    if isinstance(rho, float): rho = np.array([rho])
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - given) + x*R2*given
    w = q1*q2
    mask = (w > 0.0) & (rho == 0.0)
    w[mask] = pow(q1[mask], alpha[mask])*pow(q2[mask], 1.0 - alpha[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = alpha[mask]*pow(q1[mask], rho[mask]) + (1.0 - alpha[mask])*pow(q2[mask], rho[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = pow(w[mask], 1.0/rho[mask])
    mask = (rho > 0.0)
    w[mask] = pow(alpha[mask]*pow(q1[mask], rho[mask]) + (1.0 - alpha[mask])*pow(q2[mask], rho[mask]), 1.0/rho[mask])
    return w

b = a2max/a1max
givens = np.linspace(maxgiven, mingiven, num=num)
givens[0] = 0.999999
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
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

for alpha in alphas:

    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(18, 18))
    fig.delaxes(axs[0, 1])
    fig.delaxes(axs[0, 2])
    fig.supylabel("Partner's share of $\it{A}$", x=0.04, y=0.520, fontsize=fslabel)
    fig.supxlabel("Substitutability of $\it{A}$", x=0.525, y=0.05, fontsize=fslabel)

    if movie: fig.text(0.93, 0.02, f'alpha = {alpha}', fontsize=fstick, color='grey', ha='right')

    Q0 = Rq*pow(T0*(1.0 - alpha)/alpha, 1.0/(RR - 1.0))
    a20ss = a2max/(1.0 + Q0*b)
    Q = Rq*pow(TT*(1.0 - alpha)/alpha, 1.0/(RR - 1.0))
    a2eqss = a2max/(1.0 + Q*b)

    A = np.full([npoints, npoints], alpha)
    #A = alpha
    Zss = np.empty((0, npoints*num))
    for given in givens:
        G = np.full([npoints, npoints], given)
        Zs = np.empty((npoints, 0))
        for rho in rhos:
            Rh = np.full([npoints, npoints], rho)
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
            Z[(T > R) & (P > S)] = 0.1
            Z[(T >= R) & (P <= S)] = 0.5
            Z[(T < R) & (P < S)] = 0.9
            Z[(T == R) & (P == R)] = 1.0
            #Z = np.tril(Z, k=-1)
            Z = np.ma.masked_where(Z == 0.0, Z)
            Zs = np.append(Zs, Z, axis=1)
        Zss = np.append(Zss, Zs, axis=0)

    AA = np.full([num, num], alpha)
    Mss = [[a20ss, a20ss*R2*GG, fitness(a20ss, a20ss, GG, AA, RR)], [a2eqss, a2eqss*R2*GG, fitness(a2eqss, a2eqss, GG, AA, RR)]]

    minx = round(log_ess[0])
    maxx = round(log_ess[-1])
    xticklabels = [minx, round((minx + maxx)/2), maxx]
    yticklabels = [0.0, 0.5, 1.0]

    for axrow, letterrow in zip(axs, letters):
        for ax, letter, traitlabel in zip(axrow, letterrow, traitlabels):
            if ax.get_subplotspec().is_first_row():
                ax.set(xticks=[npoints*num, npoints*num/2, 0], yticks=[npoints*num, npoints*num/2, 0], xticklabels=[], yticklabels=yticklabels, xlim=(0, npoints*num), ylim=(npoints*num, 0))
                #ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)
                ax.text(0, 1.2, letter, fontsize=fslabel, weight='bold')
            else:
                ax.set(xticks=[-0.5, num/2, num-0.5], yticks=[num-0.5, num/2, -0.5], xticklabels=[], yticklabels=[], xlim=(-0.5, num-0.5), ylim=(num-0.5, -0.5))
                ax.text(0, 1.3, letter, fontsize=fslabel, weight='bold')
            if ax.get_subplotspec().is_last_row():
                ax.set_xticklabels(xticklabels, fontsize=fstick)
            if ax.get_subplotspec().is_first_col():
                ax.set_yticklabels(yticklabels, fontsize=fstick) 

    axs[0, 0].imshow(Zss, origin='lower', cmap='magma', vmin=0, vmax=1)

    for row, Ms in zip(axs[1:], Mss):
        for ax, M, traitvmax in zip(row, Ms, traitvmaxs):
            ax.imshow(M, cmap='magma', vmin=0, vmax=traitvmax)

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
