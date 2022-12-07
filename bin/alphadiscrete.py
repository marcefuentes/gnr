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

letters = [['a', 'b', 'c'],
            ['d', 'e', 'f'],
            ['g', 'h', 'i']]
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
vmax = 1.5

num = 21    # Number of subplot rows and columns
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

# Theory

b = a2max/a1max
givens = np.linspace(maxgiven, mingiven, num=num)
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)

# Figure properties

fslabel=24 # Label font size
fstick=16 # Tick font size

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

yellow = [1.0, 0.7, 0.1, 1.0]
red = [1.0, 0.2, 0.2, 1.0]
blue = [0.2, 0.0, 0.4, 1.0]

fig = plt.figure(figsize=(12, 12)) 
fig.supxlabel('Substitutability of $\it{A}$', x=0.52, y=0.06, fontsize=fslabel*1.35)
fig.supylabel('Partner\'s share of $\it{A}$', x=0.06, y=0.53, fontsize=fslabel*1.35, ha='center')

outer_grid = fig.add_gridspec(nrows=len(alphas), ncols=3, left=0.15, right=0.9, top=0.9, bottom=0.15)

for g, alpha, alphafolder, letterrow in zip(range(3), alphas, alphafolders, letters):
    grid = outer_grid[g, 0].subgridspec(num, num, wspace=0, hspace=0)
    axs = grid.subplots()

    # Theory

    givens[0] = 1.0
    a2 = np.linspace(0.0, a2max, num=3)
    xaxis = [1, 2, 3, 4]

    for row, given in zip(axs, givens):
        for ax, rho in zip(row, rhos):
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
                rgb = yellow
            elif (T >= R) & (P <= S) & (R - S - T + P != 0.0):
                x = (P - S)/(R - S - T + P)
                rgb = red
            elif (T > R) & (P > S):
                x = 0.0
                rgb = blue
            else:
                x = 0.0
                rgb = 'cyan'
            yaxis = [T, R, P, S]
            ax.plot(xaxis, yaxis, color=rgb, marker='o', markerfacecolor='white', linewidth=1, markersize=2)
            ax.set(xticks=[], yticks=[], xlim=(0, 5), ylim=(0.0, 2.0))
        for ax, given in zip(axs[::every, 0], givens[::every]):
            ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)
    for ax, given in zip(axs[::every, 0], givens[::every]):
        ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)
    if g == 2:
        for ax, log_es in zip(axs[-1, ::every], log_ess[::every]):
            ax.set_xlabel(round(log_es), fontsize=fstick)

    axs[0, 0].set_title(letterrow[0], fontsize=fslabel, weight='bold')

    # Simulations

    for h, folder, trait, traitlabel, traitvmax, letter in zip(range(2), folders, traits, traitlabels, traitvmaxs, letterrow[1:]):
        grid = outer_grid[g, h+1].subgridspec(1, 1)
        ax = grid.subplots()

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
        ax.set(xticks=[], yticks=[], xticklabels=[], yticklabels=[])
        if g == 2:
            xticklabels = [minx, round((minx + maxx)/2), maxx]
            ax.set_xticks([0, nc/2, nc])
            ax.set_xticklabels(xticklabels, fontsize=fstick)
            ax.xaxis.set_tick_params(length=0, width=0)
        ax.text(0, nr*1.035, letter, fontsize=fslabel, weight='bold')
        ax.imshow(df_piv, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)

plt.savefig('alphadiscrete.png', transparent=False)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
