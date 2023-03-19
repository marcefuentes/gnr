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

def figdata(given, budget, icurve):
    a2private = my.a2eq(given, AA, RR)
    w = my.fitness(a2private, a2private, given, AA, RR)
    q2_partner = a2private*my.R2
    budget_own = budget0*(1.0 - given)

    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            budgety = budget_own + q2_partner[i, j]*given
            budget[i, j].set_ydata(budgety)
            icy = my.indifference(icx, w[i, j], alpha, rho)
            icurve[i, j].set_ydata(icy)
            icurve[i, j].set_color(cm.viridis(w[i, j]/traitvmax))
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

# Figure properties

xlim=[0.0, my.a1max*my.R1]
ylim=[0.0, my.a2max*my.R2]
step = int(num/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
traitvmax = my.fitness(np.array([my.a2max]),
                       np.array([my.a2max]),
                       np.array([0.0]),
                       np.array([0.9]),
                       np.array([5.0]))

width = plotsize*len(givens)
height = plotsize
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.502,
              y=0.0,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.07,
              y=0.5,fontsize=biglabels)

outergrid = fig.add_gridspec(nrows=1,
                             ncols=len(givens),
                             left=0.15,
                             right=0.85,
                             top=0.8,
                             bottom=0.2)

for g, given in enumerate(givens):
    grid = outergrid[g].subgridspec(nrows=num,
                                    ncols=num,
                                    wspace=0,
                                    hspace=0)
    axs = grid.subplots()
    letter = ord('a') + g
    axs[0, 0].set_title(chr(letter),
                        fontsize=plotsize*5,
                        weight='bold',
                        loc='left')

    for ax in fig.get_axes():
        ax.set(xticks=[], yticks=[])
        ax.set(xlim=xlim, ylim=ylim)
    if g == 0:
        for i in range(0, num, step):
            axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                                 rotation='horizontal',
                                 horizontalalignment='right',
                                 verticalalignment='center',
                                 fontsize=ticklabels)
    for j in range(0, num, step):
        axs[-1, j].set_xlabel(f'{logess[j]:2.0f}',
                              x=0.45,
                              fontsize=ticklabels)

    budget = np.empty(axs.shape, dtype=object)
    icurve = np.empty(axs.shape, dtype=object)
    dummy_budgety = np.zeros_like(budgetx)
    dummy_icy = np.zeros_like(icx)

    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            for k in range(n_ic): 
                axs[i, j].plot(icx, ics[i, j, k], c='0.850')
            budget[i, j], = axs[i, j].plot(budgetx,
                                           dummy_budgety,
                                           c='black',
                                           alpha=0.8)
            icurve[i, j], = axs[i, j].plot(icx,
                                           dummy_icy,
                                           linewidth=4,
                                           alpha=0.8)

    axs[0, int(num/2)].set_title(f'{given*100:2.0f}%',
                                 pad=plotsize*5,
                                 fontsize=plotsize*5)
# Save figure

figdata(given, budget, icurve,)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
