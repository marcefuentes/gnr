#! /usr/bin/env python

import os
import time

import matplotlib.pyplot as plt
import numpy as np

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

givens = [1.0, 0.95, 0.5]

num = 21    # Number of subplot rows and columns
rows = givens
plotsize = 8

alphas = np.linspace(my.alphamax, my.alphamin, num=num)
logess = np.linspace(my.logesmin, my.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
high = my.a2eq(0.0, AA, RR)

xlim=[0, 5]
ylim=[0.0, 2.0]
step = int(num/2)
xaxis = [1, 2, 3, 4]
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
width = plotsize
height = plotsize*len(rows)
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(width*1.13, height))
outergrid = fig.add_gridspec(nrows=len(rows),
                             ncols=1,
                             left=0.25,
                             right=0.85)

axs = np.empty((len(rows),
                len(alphas),
                len(rhos)),
                dtype=object)
for g, row in enumerate(rows):
    grid = outergrid[g].subgridspec(nrows=num,
                                    ncols=num,
                                    wspace=0,
                                    hspace=0)
    axs[g] = grid.subplots()
    axs[g, 0, 0].set_title(chr(letter),
                              fontsize=plotsize*5,
                              weight='bold',
                              loc='left')
    letter += 1
    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            axs[g, i, j].set(xticks=[], yticks=[])
            axs[g, i, j].set(xlim=xlim, ylim=ylim)
            for axis in ['top','bottom','left','right']:
                axs[g, i, j].spines[axis].set_linewidth(0.1)
    for i in range(0, num, step):
        axs[g, i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                             rotation='horizontal',
                             horizontalalignment='right',
                             verticalalignment='center',
                             fontsize=ticklabels)
    if g == 2:
        for j in range(0, num, step):
            axs[g, -1, j].set_xlabel(f'{logess[j]:2.0f}',
                                        x=0.3,
                                        fontsize=ticklabels)

left_x = axs[0, 0, 0].get_position().x0
right_x = axs[0, -1, -1].get_position().x1
center_x = (left_x + right_x) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=0.04,
              fontsize=biglabels)

top_y = axs[0, 0, 0].get_position().y1
bottom_y = axs[-1, -1, 0].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supylabel(ylabel,
              x=0.04,
              y=center_y,
              fontsize=biglabels)

# Create lines

for g, given in enumerate(givens):

    low = my.a2eq(given, AA, RR)

    T = my.fitness(high, low, given, AA, RR)
    R = my.fitness(high, high, given, AA, RR)
    P = my.fitness(low, low, given, AA, RR)
    S = my.fitness(low, high, given, AA, RR)
    Z = my.gamecolors(T, R, P, S)
    greys = np.full(Z.shape, [0.8, 0.8, 0.8, 1.0])
    m = (Z == my.colormap['white'])
    Z[m] = greys[m]

    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            y = [T[i, j], R[i, j], P[i, j], S[i, j]]
            for line in axs[g, i, j].get_lines():
                line.remove()
            axs[g, i, j].plot(xaxis,
                           y,
                           c=Z[i, j],
                           linewidth=3,
                           marker='o',
                           markerfacecolor='white',
                           markersize=plotsize/3)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
