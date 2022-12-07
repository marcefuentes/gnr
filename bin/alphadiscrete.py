#! /usr/bin/env python

from glob import glob
from math import log
import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
import time

start_time = time.perf_counter ()

movie = False

letters = [['a', 'b', 'c', 'd'],
            ['e', 'f', 'g', 'h'],
            ['i', 'j', 'k', 'l']]
traits = ['ChooseGrainmean', 'MimicGrainmean']
traitlabels = ['Sensitivity for\nchoosing partner', 'Sensitivity for\nmimicking partner']
traitvmaxs = [1.0, 1.0]
folders = ['p', 'r']
alphafolders = ['alpha25', 'alpha50', 'alpha75']
alphas = [0.25, 0.50, 0.75]

R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0
npoints = 128
vmax = 1.5
limmatrix = a2max

num = 11    # Number of subplot rows and columns
every = int(num/2)
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0

def fitness(x, y, alpha, given, rho):
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

# No cooperation

b = a2max/a1max
givens = np.linspace(maxgiven, mingiven, num=num)
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)

# Figure properties

fslabel=36 # Label font size
fstick=24 # Tick font size

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
frames = []

fig = plt.figure(figsize=(12, 12)) 
#fig.supxlabel('Substitutability of $\it{A}$', x=0.513, y=0.05, fontsize=fslabel*1.25)
#fig.supylabel('Partner\'s share of $\it{A}$', x=0.05, y=0.493, fontsize=fslabel*1.25, ha='center')

outer_grid = fig.add_gridspec(len(alphas), 3, left=0.15, right=0.9, top=0.9, bottom=0.15)

for g, alpha, alphafolder, letterrow in zip(range(3), alphas, alphafolders, letters):
    print(alpha)
    grid = outer_grid[g, 0].subgridspec(num, num, wspace=0, hspace=0)
    axs0 = grid.subplots()
    grid = outer_grid[g, 1].subgridspec(num, num, wspace=0, hspace=0)
    axs1 = grid.subplots()

    givens[0] = 1.0
    a2 = np.linspace(0.0, a2max, num=3)
    xaxis = [1, 2, 3, 4]

    for row0, row1, given in zip(axs0, axs1, givens):
        for ax0, ax1, rho in zip(row0, row1, rhos):
            w = []
            T = fitness(a2[1], a2[0], alpha, given, rho)
            R = fitness(a2[1], a2[1], alpha, given, rho)
            P = fitness(a2[0], a2[0], alpha, given, rho)
            S = fitness(a2[0], a2[1], alpha, given, rho)
            if R < P:
                H = T
                T = S
                S = H
                H = R
                R = P
                P = H
            Su = fitness(a2[1], a2[2], alpha, given, rho)
            if Su > R:
                T = fitness(a2[2], a2[1], alpha, given, rho)
                R = fitness(a2[2], a2[2], alpha, given, rho)
                P = fitness(a2[1], a2[1], alpha, given, rho)
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
                w.append(fitness(a2[0], a, alpha, given, rho)*(1.0 - x) + fitness(a2[1], a, alpha, given, rho)*x)   
            yaxis = [T, R, P, S]
            ax0.plot(xaxis, yaxis, color=rgb, marker='o', markerfacecolor='white', linewidth=2, markersize=3)
            ax0.set(xticks=[], yticks=[], xlim=(0, 5), ylim=(0.0, 2.0))
            ax0.set_facecolor('0.200')
            ax0.set_box_aspect(1)
            ax1.plot(a2, w, linewidth=2, c=cm.magma(weq/vmax))
            ax1.set(xticks=[], yticks=[], xlim=(0.0, a2max), ylim=(0.0, 2.0))
            ax1.set_facecolor('0.200')
            ax1.set_box_aspect(1)

    axs0[0, 0].set_title('b', fontsize=fslabel, weight='bold')
    axs1[0, 0].set_title('d', fontsize=fslabel, weight='bold')
    for ax1, log_es in zip(axs1[-1, ::every], log_ess[::every]):
        ax1.set_xlabel(round(log_es), fontsize=fstick)

    grid = outer_grid[g, 2].subgridspec(1, 2)
    axs = grid.subplots()

    for ax, folder, trait, traitlabel, traitvmax, letter in zip(axs, folders, traits, traitlabels, traitvmaxs, letterrow):
        df = pd.concat(map(pd.read_csv, glob(os.path.join(alphafolder, 'discrete', folder, '*.csv'))), ignore_index=True)
        ts = df.Time.unique()
        t = ts[-1]
        df = df.loc[df.Time == t].copy()
        sgivens = np.sort(pd.unique(df.Given))[::-1]
        sess = np.sort(pd.unique(df.ES))
        df[trait] = 1.0 - df[trait]
        df_piv = pd.pivot_table(df, values=trait, index=['Given'], columns=['ES']).sort_index(axis=0, ascending=False)
        srhos = 1.0 - 1.0/sess
        nr = len(sgivens)
        nc = len(srhos)
        extent = 0, nr, 0, nc
        minx = round(log(sess[0], 2))
        maxx = round(log(sess[-1], 2))
        miny = sgivens[-1]
        maxy = sgivens[0]
        xticklabels = [minx, round((minx + maxx)/2), maxx]
        yticklabels = [miny, (miny + maxy)/2, maxy]

        ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
        ax.set_yticklabels(yticklabels, fontsize=fstick) 
        ax.set_xticklabels(xticklabels, fontsize=fstick)
        ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)
        ax.text(0, nr*1.035, letter, fontsize=fslabel, weight='bold')

        ax.imshow(df_piv, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)

plt.savefig('alphadiscrete.png', transparent=False)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
