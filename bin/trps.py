#! /usr/bin/env python

import os
import time

import matplotlib.pyplot as plt
import numpy as np

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

# Options

givens = [1.0, 0.95, 0.5]

num = 21    # Number of subplot rows and columns
plotsize = 8

# Add data to figure

def figdata(lines):

    for g, given in enumerate(givens)

        low = my.a2eq(given, AA, RR)

        T = my.fitness(high, low, given, AA, RR)
        R = my.fitness(high, high, given, AA, RR)
        P = my.fitness(low, low, given, AA, RR)
        S = my.fitness(low, high, given, AA, RR)
        Z = my.gamecolors(T, R, P, S)
        greys = np.full(Z.shape, [0.8, 0.8, 0.8, 1.0])
        m = (Z == my.colormap['white'])
        Z[m] = greys[m]
        for a, alpha in enumerate(alphas):
            for r, rho in enumerate(rhos):
                y = [T[a, r], R[a, r], P[a, r], S[a, r]]
                lines[g, a, r].set.y_data(y)
                lines[g, a, r].set_color(Z[a, r])

    return lines.flatten()

# Get data

alphas = np.linspace(my.alphamax, my.alphamin, num=num)
logess = np.linspace(my.logesmin, my.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
high = my.a2eq(0.0, AA, RR)

# Figure properties

width = plotsize
height = plotsize*len(rows)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
xlim=[0, 5]
ylim=[0.0, 2.0]
step = int(num/2)
xaxis = [1, 2, 3, 4]
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig = plt.figure(figsize=(width*1.13, height))
outergrid = fig.add_gridspec(nrows=len(rows),
                             ncols=1,
                             left=0.25,
                             right=0.85)

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
        ax.spines[axis].set_linewidth(0.1)

for g, given in enumerate(givens):
    letter = ord('a') + g
    axs[g, 0, 0].set_title(chr(letter),
                           fontsize=plotsize*5,
                           weight='bold',
                           loc='left')
    for a in range(0, num, step):
        axs[g, a, 0].set_ylabel(f'{alphas[a]:.1f}',
                                rotation='horizontal',
                                horizontalalignment='right',
                                verticalalignment='center',
                                fontsize=ticklabels)
    if g == 2:
        for r in range(0, num, step):
            axs[g, -1, r].set_xlabel(f'{logess[r]:.0f}',
                                     x=0.3,
                                     fontsize=ticklabels)

# Assign Line2D objects to lines

dummy_y = np.zeros_like(xaxis)

for g, given in enumerate(givens):
    for a, alpha in enumerate(alphas):
        for r, rho in enumerate(rhos):
            lines[g, a, r], = axs[g, a, r].plot(xaxis,
                                                dummy_y,
                                                linewidth=3,
                                                marker='o',
                                                markerfacecolor='white',
                                                markersize=plotsize/3)

# Add data and save figure

figdata(lines,)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
