#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule as my
import numpy as np
import os
import time

start_time = time.perf_counter()
this_file = os.path.basename(__file__)
file_name = this_file.split(".")[0]

givens = [1.0, 0.95, 0.5]

num = 21    # Number of subplot rows and columns
ext = 256
rows = givens
plotsize = 8

alphas = np.linspace(my.alphamax, my.alphamin, num=num)
logess = np.linspace(my.logesmin, my.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)

step = int(num/2)
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
letter = ord("a")
extent = 0, ext, 0, ext
width = plotsize
height = plotsize*len(rows)
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

fig = plt.figure(figsize=(width*1.13, height))
fig.supxlabel(xlabel,
              x=0.55,
              y=0.04,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.04,
              y=0.495,
              fontsize=biglabels)

outergrid = fig.add_gridspec(nrows=len(rows),
                             ncols=1,
                             left=0.25,
                             right=0.85)

axss = []
for g, row in enumerate(rows):
    grid = outergrid[g].subgridspec(nrows=num,
                                    ncols=num,
                                    wspace=0,
                                    hspace=0)
    axs = grid.subplots()
    axs[0, 0].set_title(chr(letter),
                        fontsize=plotsize*5,
                        weight="bold",
                        loc="left")
    letter += 1
    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            axs[i, j].set(xticks=[], yticks=[])
            for axis in ["top","bottom","left","right"]:
                axs[i, j].spines[axis].set_linewidth(0.1)
    for i in range(0, num, step):
        axs[i, 0].set_ylabel(f"{alphas[i]:3.1f}",
                             rotation="horizontal",
                             horizontalalignment="right",
                             verticalalignment="center",
                             fontsize=ticklabels)
    if g == 2:
        for j in range(0, num, step):
            axs[-1, j].set_xlabel(f"{logess[j]:2.0f}",
                                  x=0.3,
                                  fontsize=ticklabels)
    axss.append(axs)

for g, given in enumerate(givens):

    for i, alpha in enumerate(alphas):
        AA = np.full([ext, ext], alpha)
        for j, rho in enumerate(rhos):

            xmin = 0.0
            xmax = 1.0
            ymin = 0.0
            ymax = 1.0
            x = np.linspace(xmin, xmax, num=ext)
            y = np.linspace(ymax, ymin, num=ext)
            X, Y = np.meshgrid(x, y)
            RR = np.full([ext, ext], rho)
            T = my.fitness(Y, X, given, AA, RR)
            R = my.fitness(Y, Y, given, AA, RR)
            P = my.fitness(X, X, given, AA, RR)
            S = my.fitness(X, Y, given, AA, RR)

            Z = my.gamecolors(T, R, P, S)
            Z[X >= Y] = [0.9, 0.9, 0.9, 1.0]
            axss[g][i][j].imshow(Z, extent=extent)

plt.savefig(file_name + ".png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
