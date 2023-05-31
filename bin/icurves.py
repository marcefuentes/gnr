#! /usr/bin/env python

import os
import time

from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

import mymodule as my

start_time = time.perf_counter()
this_file = os.path.basename(__file__)
file_name = this_file.split('.')[0]

# Options

givens = [0.0, 0.5, 0.95]

vmax = my.wmax
num = 3     # Number of subplot rows & columns
numaB = 256 # Number of points along each curve
n_ic = 5    # Number of indifference curves

plotsize = 6

# Add data to figure

def init(budgets, icurves):

    for g, given in enumerate(givens):
        aBprivate = my.aBeq(given, AA, RR)
        w = my.fitness(aBprivate, aBprivate, given, AA, RR)
        qB_partner = aBprivate*my.RB
        budget_own = budget0*(1.0 - given)

        for a, alpha in enumerate(alphas):
            for r, rho in enumerate(rhos):
                budgety = budget_own + qB_partner[a, r]*given
                budgets[g, a, r].set_ydata(budgety)
                icy = my.indifference(icx, w[a, r], alpha, rho)
                icurves[g, a, r].set_ydata(icy)
                color = cm.viridis(w[a, r]/vmax)
                icurves[g, a, r].set_color(color)

    return np.concatenate([budgets.flatten(), icurves.flatten()])

# Data

alphas = np.linspace(my.alphamax, my.alphamin, num=num)
logess = np.linspace(my.logesmin, my.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
a1 = np.array([0.0, my.aAmax])
budgetx = a1*my.RA
budget0 = (my.aBmax - my.b*a1)*my.RB
icx = np.linspace(0.001*my.RA,
                  (my.aAmax - 0.001)*my.RA,
                  num=numaB)
RR, AA = np.meshgrid(rhos, alphas)
ws = np.linspace(2.0/(n_ic + 1), 2.0*n_ic/(n_ic + 1), num=n_ic)
ics = np.zeros((num, num, n_ic, numaB))
for i, alpha in enumerate(alphas):
    for j, rho in enumerate(rhos):
        for k, w in enumerate(ws):
            ics[i, j, k] = my.indifference(icx, w, alpha, rho)

# Figure properties

width = plotsize*len(givens)
height = plotsize
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Influence of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
xlim=[0.0, my.aAmax*my.RA]
ylim=[0.0, my.aBmax*my.RB] 
step = int(num/2)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=1,
                             ncols=len(givens),
                             left=0.15,
                             right=0.85,
                             top=0.8,
                             bottom=0.2)
axs = np.empty((len(givens),
                len(alphas),
                len(rhos)),
               dtype=object)

for g, given in enumerate(givens):
    grid = outergrid[g].subgridspec(nrows=num,
                                    ncols=num,
                                    wspace=0,
                                    hspace=0)
    axs[g] = grid.subplots()

left_x = axs[0, 0, 0].get_position().x0
right_x = axs[-1, -1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0, 0].get_position().y1
bottom_y = axs[-1, -1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.2,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=left_x*0.5,
              y=center_y,
              fontsize=biglabels)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    ax.set(xlim=xlim, ylim=ylim)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.2)

for g, given in enumerate(givens):
    letter = ord('a') + g
    axs[g, 0, 0].set_title(chr(letter),
                           fontsize=plotsize*5,
                           weight='bold',
                           loc='left')
    axs[g, 0, int(num/2)].set_title(f'{given*100:.0f}%',
                                    pad=plotsize*5,
                                    fontsize=plotsize*5)
    if g == 0:
        for a in range(0, num, step):
            axs[g, a, 0].set_ylabel(f'{alphas[a]:.1f}',
                                    rotation='horizontal',
                                    horizontalalignment='right',
                                    verticalalignment='center',
                                    fontsize=ticklabels)
    for r in range(0, num, step):
        axs[g, -1, r].set_xlabel(f'{logess[r]:.0f}',
                                 fontsize=ticklabels)

# Assign axs objects to variables
# (Line2D)

budgets = np.empty_like(axs)
icurves = np.empty_like(axs)
dummy_budgety = np.zeros_like(budgetx)
dummy_icy = np.zeros_like(icx)

for g, given in enumerate(givens):
    for a, alpha in enumerate(alphas):
        for r, rho in enumerate(rhos):
            for c in range(n_ic): 
                axs[g, a, r].plot(icx, ics[a, r, c], c='0.850')
            budgets[g, a, r], = axs[g, a, r].plot(budgetx,
                                                  dummy_budgety,
                                                  c='black',
                                                  alpha=0.8)
            icurves[g, a, r], = axs[g, a, r].plot(icx,
                                                  dummy_icy,
                                                  linewidth=4,
                                                  alpha=0.8)

# Add data and save figure

init(budgets, icurves,)

plt.savefig(file_name + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
