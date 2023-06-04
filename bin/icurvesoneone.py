#! /usr/bin/env python

import os
import time

from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

import mymodule as my

start_time = time.time()
this_file = os.path.basename(__file__)
file_name = os.path.splitext(this_file)[0]

# Options

numaB = 256 # Number of points along each curve
n_ic = 10   # Number of indifference curves
plotsize = 8.0 # Size of figure
xlabel = 'Quantity of A'
ylabel = 'Quantity of B'

# Data

alpha = 0.7
loges = -1.0
rho = 1.0 - 1.0/pow(2.0, loges)
ws = np.linspace(2.0/(n_ic + 1), 2.0*n_ic/(n_ic + 1), num=n_ic)
x = np.linspace(0.001*my.RA, (my.aAmax - 0.001)*my.RA, num=numaB)
ys = np.zeros((n_ic, numaB))
for i, w in enumerate(ws):
    ys[i] = my.indifference(x, w, alpha, rho)

# Figure properties

xticks = [0.0, 1.0, 2.0]
yticks = [0.0, 1.0, 2.0]
xticklabels = [0.0, 0.5, 1.0]
yticklabels = [0.0, 0.5, 1.0]
width = plotsize
height = plotsize
biglabels = plotsize * 4
ticklabels = plotsize * 3

# Create figure

fig = plt.figure(figsize=(width, height))
ax = fig.add_subplot(111)
ax.set_xlabel(xlabel, fontsize=biglabels)
ax.set_ylabel(ylabel, fontsize=biglabels)
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xticklabels(xticklabels, fontsize=ticklabels)
ax.set_yticklabels(yticklabels, fontsize=ticklabels)
ax.set_xlim(0.0, 2.0)
ax.set_ylim(0.0, 2.0)
# Map of indifference curves
for i in range(n_ic):
    ax.plot(x, ys[i], linewidth=2.0, color=cm.viridis(ws[i]/1.5), alpha=0.2)
# Budget line
ax.plot([2, 0], [0, 2], linewidth=2.0, color='black', alpha=0.5)
# Equilibrium indifference curve
aBprivate = my.aBeq(0.0, alpha, rho)
w = my.fitness(aBprivate, aBprivate, 0.0, alpha, rho)
y = my.indifference(x, w, alpha, rho)
ax.plot(x, y, linewidth=2.0, color=cm.viridis(w/1.5), alpha=1.0)

# Save figure

fig.savefig(file_name + '.png', bbox_inches='tight')
plt.close(fig)
end_time = time.time()
elapsed_time = end_time - start_time
print(f"\nTime elapsed: {elapsed_time:.2f} seconds")

