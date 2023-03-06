#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

given = 0.95
num = 5    # Number of subplot rows & columns
ext = 256
distances = [0.2, 0.5, 0.8]

plotsize = 9

alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=num)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
x = np.linspace(0.0, mymodule,a2max, num=ext)
y = np.flip(x)
X, Y = np.meshgrid(x, y)
maskxy = (X >= Y)

step = int(num/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
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
    #ax.set(xlim=xlim, ylim=ylim)
for i in range(0, num, step):
    axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                         rotation='horizontal',
                         horizontalalignment='right',
                         verticalalignment='center',
                         fontsize=ticklabels)
    axs[-1, i].set_xlabel(f'{logess[i]:2.0f}',
                          fontsize=ticklabels)

for i, alpha in enumerate(alphas):
    AA = np.full(3, alpha)
    AAA = np.full([ext, ext], alpha)
    for j, rho in enumerate(rhos):
        RR = np.full(3, rho)
        a2private = mymodule.a2eq(given, AA, RR)
        a2zero = mymodule.a2eq(0.0, AA, RR)
        x = distances*a2private
        y = a2zero + distances*(mymodule.a2max - a2zero)
        
        T = mymodule.fitness(y, x, given, AA, RR)
        R = mymodule.fitness(y, y, given, AA, RR)
        P = mymodule.fitness(x, x, given, AA, RR)
        S = mymodule.fitness(x, y, given, AA, RR)

        Z = mymodule.gamecolors(T, R, P, S)
        axs[i][j].scatter(x*ext, y*ext,
                          marker='o',
                          s=markersize,
                          color=Z)
        Z = mymodule.nodilemmacolorsg(T, R, P, S)
        axs[i][j].scatter(x*ext, y*ext,
                          marker='o',
                          s=markersize,
                          color=Z)

        RRR = np.full([ext, ext], rho)
        T = mymodule.fitness(Y, X, given, AAA, RRR)
        R = mymodule.fitness(Y, Y, given, AAA, RRR)
        P = mymodule.fitness(X, X, given, AAA, RRR)
        S = mymodule.fitness(X, Y, given, AAA, RRR)

        Z = mymodule.gamecolors(T, R, P, S)
        Z[maskxy] = [0.9, 0.9, 0.9, 1.0]
        axs[i][j].imshow(Z, extent=extent, alpha=0.2)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
