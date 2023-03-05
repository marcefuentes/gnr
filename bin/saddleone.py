#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

given = 0.999
alpha = 0.5
loges = -5.0
ext = 1024
distances = np.array([0.2, 0.5, 0.8])

plotsize = 6

rho = 1.0 - 1.0/pow(2, loges)
xmin = 0.0
xmax = mymodule.a2eq(given, alpha, rho)
ymin = mymodule.a2eq(0.0, alpha, rho)
ymax = mymodule.a2max
x = np.linspace(xmin, xmax, num=ext)
y = np.linspace(ymax, ymin, num=ext)
X, Y = np.meshgrid(x, y)

xlabel = 'Effort to get $\it{B}$'
ylabel = 'Effort to get $\it{B}$'
extent = 0, ext, 0, ext
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

a2eq = mymodule.a2eq(given, alpha, rho)
a2social = mymodule.a2eq(0.0, alpha, rho)
x = distances*a2eq
y = a2social + distances*(mymodule.a2max - a2social)

T = mymodule.fitness(y, x, given, alpha, rho)
R = mymodule.fitness(y, y, given, alpha, rho)
P = mymodule.fitness(x, x, given, alpha, rho)
S = mymodule.fitness(x, y, given, alpha, rho)

Z = mymodule.gamecolors(T, R, P, S)
ax.scatter(distances*ext, distances*ext,
           marker='o',
           s=markersize,
           color='white')
Z = mymodule.nodilemmacolorsg(T, R, P, S)
ax.scatter(distances*ext, distances*ext,
           marker='o',
           s=markersize,
           color='white')

T = mymodule.fitness(Y, X, given, alpha, rho)
R = mymodule.fitness(Y, Y, given, alpha, rho)
P = mymodule.fitness(X, X, given, alpha, rho)
S = mymodule.fitness(X, Y, given, alpha, rho)

Z = mymodule.gamecolors(T, R, P, S)
ax.imshow(Z, extent=extent, alpha=0.99)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
