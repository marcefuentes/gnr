#! /usr/bin/env python

import os
import time

from matplotlib import cm
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import numpy as np

import mymodule as my

start_time = time.perf_counter()
this_file = os.path.basename(__file__)
file_name = this_file.split(".")[0]

# Options

givens = np.linspace(0.0, 1.0, num=21)
givens = [0.0]

num = 3     # Number of subplot rows & columns
numa2 = 256 # Number of points along each curve
n_ic = 5    # Number of indifference curves

plotsize = 6

# Add data to figure

def update(given, budgets, icurves):
    a2private = my.a2eq(given, AA, RR)
    w = my.fitness(a2private, a2private, given, AA, RR)
    qB_partner = a2private*my.RB
    budget_own = budget0*(1.0 - given)

    for a, alpha in enumerate(alphas):
        for r, rho in enumerate(rhos):
            budgety = budget_own + qB_partner[a, r]*given
            budgets[0, a, r].set_ydata(budgety)
            icy = my.indifference(icx, w[a, r], alpha, rho)
            icurves[0, a, r].set_ydata(icy)
            color = cm.Reds(norm(w[a, r]/my.wmax))
            icurves[0, a, r].set_color(color)
            a2p = a2private[a, r]/2.0
            landscape = my.fitness(a2p, icx/2.0, given, alpha, rho)
            icurves[1, a, r].set_ydata(landscape)

    return np.concatenate([budgets.flatten(), icurves.flatten()])

# Data

alphas = np.linspace(my.alphamax, my.alphamin, num=num)
logess = np.linspace(my.logesmin, my.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
a1 = np.array([0.0, my.aAmax])
budgetx = a1*my.RA
budget0 = (my.a2max - my.b*a1)*my.RB
icx = np.linspace(0.001*my.RA,
                  (my.aAmax - 0.001)*my.RA,
                  num=numa2)
RR, AA = np.meshgrid(rhos, alphas)
ws = np.linspace(2.0/(n_ic + 1), 2.0*n_ic/(n_ic + 1), num=n_ic)
ics = np.zeros((num, num, n_ic, numa2))
for i, alpha in enumerate(alphas):
    for j, rho in enumerate(rhos):
        for k, w in enumerate(ws):
            ics[i, j, k] = my.indifference(icx, w, alpha, rho)
norm = Normalize(vmin=0, vmax=1)

# Figure properties

width = plotsize*2
height = plotsize
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
biglabel = plotsize*4
ticklabel = plotsize*3
xlim=[0.0, my.aAmax*my.RA]
ylim=[0.0, my.a2max*my.RB]
step = int(num/2)
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Create figure

axs = np.empty((2,
                len(alphas),
                len(rhos)),
               dtype=object)

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=1,
                             ncols=2,
                             left=0.15,
                             right=0.85,
                             top=0.8,
                             bottom=0.2)
for g in range(2):
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
              y=bottom_y*0.3,
              fontsize=biglabel)
fig.supylabel(ylabel,
              x=left_x*0.4,
              y=center_y,
              fontsize=biglabel)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    ax.set(xlim=xlim, ylim=ylim)
    for axis in ["top","bottom","left","right"]:
        ax.spines[axis].set_linewidth(0.2)

for g in range(2):
    letter = ord("a") + g
    axs[g, 0, 0].set_title(chr(letter),
                           fontsize=plotsize*5,
                           weight="bold",
                           loc="left")
    if g == 0:
        for a in range(0, num, step):
            axs[g, a, 0].set_ylabel(f"{alphas[a]:.1f}",
                                    rotation="horizontal",
                                    horizontalalignment="right",
                                    verticalalignment="center",
                                    fontsize=ticklabel)
    for r in range(0, num, step):
        axs[g, -1, r].set_xlabel(f"{logess[r]:.0f}",
                                 fontsize=ticklabel)

# Assign axs objects to variables
# (Line2D)

budgets = np.empty_like(axs)
icurves = np.empty_like(axs)
dummy_budgety = np.full_like(budgetx, -1.0)
dummy_icy = np.zeros_like(icx)

for g in range(2):
    for a, alpha in enumerate(alphas):
        for r, rho in enumerate(rhos):
            if g == 0:
                for c in range(n_ic): 
                    axs[0, a, r].plot(icx, ics[a, r, c], c="0.850")
            budgets[g, a, r], = axs[g, a, r].plot(budgetx,
                                                  dummy_budgety,
                                                  c="0.300",
                                                  linewidth=4,
                                                  alpha=0.8)
            icurves[g, a, r], = axs[g, a, r].plot(icx,
                                                  dummy_icy,
                                                  linewidth=4,
                                                  alpha=0.8)

# Add colorbar
axins = inset_axes(axs[0, -1, -1],
                   width="5%",
                   height="100%",
                   loc="upper right",
                   bbox_to_anchor=(880, 200, 200, 200),
                   borderpad=0)
cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap="Reds"),
                    cax=axins,
                    ticks=[0, 0.5, 1])
cbar.ax.tick_params(labelsize=ticklabel)
cbar.outline.set_linewidth(0.2)

# Add data and save figure

if len(givens) > 1:
    ani = FuncAnimation(fig,
                        update,
                        frames=givens,
                        fargs=(budgets, icurves,),
                        blit=True)
    ani.save(f"{file_name}.mp4", writer="ffmpeg", fps=10)
else:
    update(givens[0], budgets, icurves,)
    plt.savefig(f"{file_name}.png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
