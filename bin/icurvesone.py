#! /usr/bin/env python

import os
import time

from matplotlib import cm
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

# Options

givens = np.linspace(0.0, 1.0, num=21)
#givens = [0.95]

vmax = 2.0
num = 3     # Number of subplot rows & columns
numa2 = 256 # Number of points along each curve
n_ic = 5    # Number of indifference curves

plotsize = 6

# Add data to figure

def update(given, budgets, icurves):
    a2private = my.a2eq(given, AA, RR)
    w = my.fitness(a2private, a2private, given, AA, RR)
    q2_partner = a2private*my.R2
    budget_own = budget0*(1.0 - given)
    axs[0, int(num/2)].title.set_text(f'{given*100:.0f}%')

    for a, alpha in enumerate(alphas):
        for r, rho in enumerate(rhos):
            budgety = budget_own + q2_partner[a, r]*given
            budgets[a, r].set_ydata(budgety)
            icy = my.indifference(icx, w[a, r], alpha, rho)
            icurves[a, r].set_ydata(icy)
            color = cm.viridis(w[a, r]/vmax)
            icurves[a, r].set_color(color)

    return np.concatenate([budgets.flatten(), icurves.flatten()])

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

width = plotsize
height = plotsize
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
xlim=[0.0, my.a1max*my.R1]
ylim=[0.0, my.a2max*my.R2]
step = int(num/2)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig = plt.figure(figsize=(width, height))
grid = fig.add_gridspec(nrows=num,
                        ncols=num,
                        left=0.22,
                        right=0.9,
                        top=0.86,
                        bottom=0.176,
                        wspace=0,
                        hspace=0)

axs = grid.subplots()

left_x = axs[0, 0].get_position().x0
right_x = axs[-1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0].get_position().y1
bottom_y = axs[-1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.2,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=left_x*0.2,
              y=center_y,
              fontsize=biglabels)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    ax.set(xlim=xlim, ylim=ylim)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.2)
axs[0, int(num/2)].set_title(f'0',
                             pad=plotsize*5,
                             fontsize=plotsize*5)
for a in range(0, num, step):
    axs[a, 0].set_ylabel(f'{alphas[a]:.1f}',
                         rotation='horizontal',
                         horizontalalignment='right',
                         verticalalignment='center',
                         fontsize=ticklabels)
for r in range(0, num, step):
    axs[-1, r].set_xlabel(f'{logess[r]:.0f}',
                          fontsize=ticklabels)

# Assign axs objects to variables
# (Line2D objects to lines)

budgets = np.empty(axs.shape, dtype=object)
icurves = np.empty(axs.shape, dtype=object)
dummy_budgety = np.zeros_like(budgetx)
dummy_icy = np.zeros_like(icx)

for a, alpha in enumerate(alphas):
    for r, rho in enumerate(rhos):
        for c in range(n_ic): 
            axs[a, r].plot(icx, ics[a, r, c], c='0.850')
        budgets[a, r], = axs[a, r].plot(budgetx,
                                        dummy_budgety,
                                        c='black',
                                        alpha=0.8)
        icurves[a, r], = axs[a, r].plot(icx,
                                        dummy_icy,
                                        linewidth=4,
                                        alpha=0.8)

# Add data and save figure

if len(givens) > 1:
    ani = FuncAnimation(fig, update, frames=givens, fargs=(budgets, icurves,), blit=True)
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    update(givens[0], budgets, icurves,)
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
