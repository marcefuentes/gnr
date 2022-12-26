#! /usr/bin/env python

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
minloges = -5.0
maxloges = 5.0
mingiven = 0.95
maxgiven = 0.95

num = 5    # Number of subplot rows & columns
npoints = 128
n_ic = 5    # Number of indifference curves
ngiven = 21
vmax = 1.8

filename = 'icurves_es_alpha'
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

every = int(num/2)

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

def indifference(q, w, alpha, rho):
    if rho == 0.0:
        q2 = np.piecewise(q, [q == 0.0, q > 0.0], [1000.0, lambda i: pow(w/pow(i, 1.0 - alpha), 1.0/alpha)])
    elif rho < 0.0:
        q2 = np.piecewise(q, [q == 0.0, q > 0.0], [1000.0, lambda i: np.piecewise(i, [pow(w, rho) <= (1.0 - alpha)*pow(i, rho), pow(w, rho) > (1.0 - alpha)*pow(i, rho)], [1000.0, lambda j: pow((pow(w, rho) - (1.0 - alpha)*pow(j, rho))/alpha, 1.0/rho)])])
    else:
        q2 = np.piecewise(q, [pow(w, rho) <= (1.0 - alpha)*pow(q, rho), pow(w, rho) > (1.0 - alpha)*pow(q, rho)], [-0.1, lambda i: pow((pow(w, rho) - (1.0 - alpha)*pow(i, rho))/alpha, 1.0/rho)])
    return q2

if mingiven != maxgiven:
    movie = True
    givens = np.linspace(mingiven, maxgiven, num=ngiven)
    frames = []
else:
    movie = False 
    givens = np.array([mingiven])

Rq = R2/R1
b = a2max/a1max
nc = num
nr = num
a1_budget = np.linspace(0.0, a1max, num=3)
q2_budget = (a2max - b*a1_budget)*R2
q1_budget = a1_budget*R1
q1_ic = np.linspace(0.0, a1max*R1, num=npoints)
if givens[-1] == 0.9999999:
    givens[-1] = 0.9999999
alphas = np.linspace(maxalpha, minalpha, num=nr)
logess = np.linspace(minloges, maxloges, num=nc)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
wis = np.linspace(2.0/(n_ic + 1), 2.0*n_ic/(n_ic + 1), num=n_ic)

xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'

icsss = []
for alpha in alphas:
    icss = []
    for rho in rhos:
        ics = []
        for w in wis:
            ics.append(indifference(q1_ic, w, alpha, rho))
        icss.append(ics)
    icsss.append(icss)

fig = plt.figure(figsize=(8, 8))
fig.supxlabel(xlabel, x=0.555, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.05, y=0.52, fontsize=fslabel)
grid = fig.add_gridspec(nrows=nr, ncols=nc, left=0.22, right=0.9, top=0.86, bottom=0.176, wspace=0, hspace=0)
axs = grid.subplots()

for row in axs:
    for ax in row:
        ax.set(xticks=[], yticks=[], xlim=(0.0, a1max*R1), ylim=(0.0, a2max*R2))
for ax, loges in zip(axs[-1, ::every], logess[::every]):
    ax.set_xlabel(round(loges), fontsize=fstick)
for ax, alpha in zip(axs[::every, 0], alphas[::every]):
    ax.set_ylabel(f'{alpha:1.1f}', rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

for given in givens:

    if movie:
        text = fig.text(0.80, 0.90, f'given: {given:4.2f}', fontsize=fstick, color='grey', ha='right')

    T = b*Rq*(1.0 - given)
    Q = Rq*pow(T*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2ss = a2max/(1.0 + Q*b)
    wss = fitness(a2ss, a2ss, given, AA, RR)
    q2ss = a2ss*R2

    for row, alpha, q2s, ws, icss in zip(axs, alphas, q2ss, wss, icsss):
        budget0 = q2_budget*(1.0 - given)
        for ax, rho, ics, q2eq, weq in zip(row, rhos, icss, q2s, ws):
            for line in ax.get_lines():
                line.remove()
            for ic in ics:
                ax.plot(q1_ic, ic, c='0.850')
            budget = budget0 + q2eq*given
            ax.plot(q1_budget, budget, c='green', alpha=0.8)
            ax.plot(q1_ic, indifference(q1_ic, weq, alpha, rho), linewidth=4, alpha= 0.8, c=cm.magma(weq/vmax))

    if movie:
        plt.savefig('temp.png', transparent=False)
        frames.append(iio.imread('temp.png'))
        os.remove('temp.png')
        text.remove()
    else:
        plt.savefig(filename + '.png', transparent=False)

plt.close()

if movie:
    iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
