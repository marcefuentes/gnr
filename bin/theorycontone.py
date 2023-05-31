#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule as my
import numpy as np
import os
import time

start_time = time.perf_counter()
this_file = os.path.basename(__file__)
filename = this_file.split('.')[0]

given = 0.95
alpha = 0.5
loges = 0.0
ext = 1025
plotsize = 6

rho = 1.0 - 1.0/pow(2, loges)
xmin = 0.0
xmax = my.aBeq(given, alpha, rho)
ymin = my.aBeq(0.0, alpha, rho)
ymax = my.aBmax
x = np.linspace(xmin, xmax, num=ext)
y = np.linspace(ymax, ymin, num=ext)
X, Y = np.meshgrid(x, y)
G = np.full([ext, ext, 4], my.colormap['transparent'])
G[X >= Y] = [0.9, 0.9, 0.9, 1.0]
T = my.fitness(Y, X, given, alpha, rho)
R = my.fitness(Y, Y, given, alpha, rho)
P = my.fitness(X, X, given, alpha, rho)
S = my.fitness(X, Y, given, alpha, rho)
Z = my.gamecolors(T, R, P, S)
Z[X >= Y] = [0.9, 0.9, 0.9, 1.0]

xlabel = 'Effort to get $\it{B}$'
ylabel = 'Effort to get $\it{B}$'
xticks = [0, ext/2.0, ext-0.5]
yticks = [0, ext/2.0, ext-0.5]
xticklabels = [f'{round(xmin):4.2f}',
               f'{round((xmax - xmin)/2.0, 3):4.2f}',
               f'{round(xmax, 3):4.2f}']
yticklabels = [f'{round(ymin, 3):4.2f}',
               f'{round((ymax + ymin)/2.0, 3):4.2f}',
               f'{round(ymax, 3):4.2f}']
markersize = plotsize*30
width = plotsize
height = plotsize
biglabels = plotsize*4
ticklabels = plotsize*3
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig, ax = plt.subplots(figsize=(width, height))
plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

ax.set_xlabel(xlabel,
              labelpad=10,
              fontsize=biglabels)
ax.set_ylabel(ylabel,
              labelpad=10,
              fontsize=biglabels)
ax.set(xticks=xticks, yticks=yticks)
ax.set_xticklabels(xticklabels, fontsize=ticklabels)
ax.set_yticklabels(yticklabels, fontsize=ticklabels)

ax.imshow(Z)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
