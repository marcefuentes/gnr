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
num = 5    # Number of subplot rows & columns
ext = 256
distances = np.array([0.2, 0.5, 0.8])

plotsize = 9

alphas = np.linspace(my.alphamax, my.alphamin, num=num)
logess = np.linspace(my.logesmin, my.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)

step = int(num/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Influence of $\it{B}$'
extent = 0, ext, 0, ext
markersize = plotsize*4
width = plotsize
height = plotsize
biglabels = plotsize*4
ticklabels = plotsize*3
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.55,
              y=0.05,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.04,
              y=0.55,
              fontsize=biglabels)

grid = fig.add_gridspec(nrows=num,
                        ncols=num,
                        wspace=0,
                        hspace=0,
                        left=0.20,
                        right=0.9,
                        top=0.9,
                        bottom=0.20)
axs = grid.subplots()
for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
for i in range(0, num, step):
    axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                         rotation='horizontal',
                         horizontalalignment='right',
                         verticalalignment='center',
                         fontsize=ticklabels)
    axs[-1, i].set_xlabel(f'{logess[i]:2.0f}',
                          fontsize=ticklabels)

for i, alpha in enumerate(alphas):
    AA = np.full([ext, ext], alpha)
    for j, rho in enumerate(rhos):
        aBprivate = my.aBeq(given, alpha, rho)
        aBsocial = my.aBeq(0.0, alpha, rho)
        x = distances*aBprivate
        y = aBsocial + distances*(my.aBmax - aBsocial)
        
        T = my.fitness(y, x, given, alpha, rho)
        R = my.fitness(y, y, given, alpha, rho)
        P = my.fitness(x, x, given, alpha, rho)
        S = my.fitness(x, y, given, alpha, rho)

        Z = my.gamecolors(T, R, P, S)
        axs[i][j].scatter(distances*ext, distances*ext,
                          marker='o',
                          s=markersize,
                          color=Z)
        Z = my.nodilemmacolorsg(T, R, P, S)
        axs[i][j].scatter(distances*ext, distances*ext,
                          marker='o',
                          s=markersize,
                          color=Z)

        x = np.linspace(0.0, aBprivate, num=ext)
        y = np.linspace(my.aBmax, aBsocial, num=ext)
        X, Y = np.meshgrid(x, y)
        RR = np.full([ext, ext], rho)
        T = my.fitness(Y, X, given, AA, RR)
        R = my.fitness(Y, Y, given, AA, RR)
        P = my.fitness(X, X, given, AA, RR)
        S = my.fitness(X, Y, given, AA, RR)

        Z = my.gamecolors(T, R, P, S)
        axs[i][j].imshow(Z, extent=extent, alpha=0.2)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
