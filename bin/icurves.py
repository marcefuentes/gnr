#! /usr/bin/env python

import os
import time

from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

# Options

givens = [0.0, 0.5, 0.95]

num = 3     # Number of subplot rows & columns
numa2 = 256 # Number of points along each curve
n_ic = 5    # Number of indifference curves

plotsize = 6

# Add data to figure

def figdata(budget, icurve):

    for i, given in enumerate(givens):
        a2private = my.a2eq(given, AA, RR)
        w = my.fitness(a2private, a2private, given, AA, RR)
        q2_partner = a2private*my.R2
        budget_own = budget0*(1.0 - given)

        for j, alpha in enumerate(alphas):
            for k, rho in enumerate(rhos):
                budgety = budget_own + q2_partner[j, k]*given
                budget[i, j, k].set_ydata(budgety)
                icy = my.indifference(icx, w[j, k], alpha, rho)
                icurve[i, j, k].set_ydata(icy)
                color = cm.viridis(w[j, k]/traitvmax)
                icurve[i, j, k].set_color(color)

    return np.concatenate([budget.flatten(), icurve.flatten()])

# Get data

alphas = np.linspace(my.alphamax, my.alphamin, num=num)
logess = np.linspace(my.logesmin, my.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
a1 = np.array([0.0, my.a1max])
budgetx = a1*my.R1
budget0 = (my.a2max - my.b*a1)*my.R2
icx = np.linspace(0.001*my.R1,
                  (my.a1max - 0.001)*my.R1,
                  num=numa2)
RR, AA = np.meshgrid(rhos, alphas)
ws = np.linspace(2.0/(n_ic + 1), 2.0*n_ic/(n_ic + 1), num=n_ic)
ics = np.empty((num, num, n_ic, numa2), dtype=np.float64)
for i, alpha in enumerate(alphas):
    for j, rho in enumerate(rhos):
        for k, w in enumerate(ws):
            ics[i, j, k] = my.indifference(icx, w, alpha, rho)
outer_columns = givens
inner_columns = logess
inner_rows = alphas

# Figure properties

width = plotsize*len(outer_columns)
height = plotsize
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
xlim=[0.0, my.a1max*my.R1]
ylim=[0.0, my.a2max*my.R2]
step = int(num/2)
traitvmax = my.fitness(np.array([my.a2max]),
                       np.array([my.a2max]),
                       np.array([0.0]),
                       np.array([0.9]),
                       np.array([5.0]))

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=1,
                             ncols=len(outer_columns),
                             left=0.15,
                             right=0.85,
                             top=0.8,
                             bottom=0.2)
axs = np.empty((len(outer_columns),
                len(inner_rows),
                len(inner_columns)),
                dtype=object)
budget = np.empty(axs.shape, dtype=object)
icurve = np.empty(axs.shape, dtype=object)

for i, outer_column in enumerate(outer_columns):
    grid = outergrid[i].subgridspec(nrows=num,
                                    ncols=num,
                                    wspace=0,
                                    hspace=0)
    axs[i] = grid.subplots()

left_x = axs[0, 0, 0].get_position().x0
right_x = axs[-1, -1, -1].get_position().x1
center_x = (left_x + right_x) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=0.0,
              fontsize=biglabels)

top_y = axs[0, 0, 0].get_position().y1
bottom_y = axs[-1, -1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supylabel(ylabel,
              x=0.07,
              y=center_y,
              fontsize=biglabels)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    ax.set(xlim=xlim, ylim=ylim)

for i, outer_column in enumerate(outer_columns):
    letter = ord('a') + i
    axs[i, 0, 0].set_title(chr(letter),
                           fontsize=plotsize*5,
                           weight='bold',
                           loc='left')
    axs[i, 0, int(num/2)].set_title(f'{outer_columns[i]*100:.0f}%',
                                    pad=plotsize*5,
                                    fontsize=plotsize*5)
    if i == 0:
        for j in range(0, num, step):
            axs[i, j, 0].set_ylabel(f'{inner_rows[j]:.1f}',
                                    rotation='horizontal',
                                    horizontalalignment='right',
                                    verticalalignment='center',
                                    fontsize=ticklabels)
    for k in range(0, num, step):
        axs[i, -1, k].set_xlabel(f'{inner_columns[k]:.0f}',
                                 x=0.45,
                                 fontsize=ticklabels)

# Assign lines to axs

dummy_budgety = np.zeros_like(budgetx)
dummy_icy = np.zeros_like(icx)

for i, given in enumerate(givens):
    for j, alpha in enumerate(alphas):
        for k, rho in enumerate(rhos):
            for l in range(n_ic): 
                axs[i, j, k].plot(icx, ics[j, k, l], c='0.850')
            budget[i, j, k], = axs[i, j, k].plot(budgetx,
                                                 dummy_budgety,
                                                 c='black',
                                                 alpha=0.8)
            icurve[i, j, k], = axs[i, j, k].plot(icx,
                                                 dummy_icy,
                                                 linewidth=4,
                                                 alpha=0.8)

# Add data and save figure

figdata(budget, icurve,)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
