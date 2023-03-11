#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule as my
import numpy as np
import os
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

givens = [1.0, 0.95, 0.5]
titles = []
for given in givens:
    titles.append(f'{given*100:2.0f}%')

num = 21    # Number of subplot rows and columns
rows = givens
plotsize = 8

alphas = np.linspace(my.alphamax, my.alphamin, num=num)
logess = np.linspace(my.logesmin, my.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
highs = [] 
for given in givens:
    highs.append(my.a2eq(0.0, AA, RR))

xlim=[0, 5]
ylim=[0.0, 2.0]
step = int(num/2)
xaxis = [1, 2, 3, 4]
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
width = plotsize
height = plotsize*len(titles)
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(width*1.25, height))
fig.supxlabel(xlabel,
              x=0.502,
              y=0.05,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.05,
              y=0.5,
              fontsize=biglabels)

outergrid = fig.add_gridspec(nrows=len(rows),
                             ncols=1,
                             left=0.20,
                             right=0.80)
                             #top=0.8,
                             #bottom=0.2)

axss = []
for g, row in enumerate(rows):
    grid = outergrid[g].subgridspec(nrows=num,
                                    ncols=num,
                                    wspace=0,
                                    hspace=0)
    axs = grid.subplots()
    #axs[int(num/2), 0].set_title(titles[g],
    #                             pad=plotsize*5,
    #                             fontsize=plotsize*5)
    axs[0, 0].set_title(chr(letter),
                        fontsize=plotsize*5,
                        weight='bold',
                        loc='left')
    letter += 1

    for ax in fig.get_axes():
        ax.set(xticks=[], yticks=[])
        ax.set(xlim=xlim, ylim=ylim)
    if g == 0:
        for i in range(0, num, step):
            axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                                 rotation='horizontal',
                                 horizontalalignment='right',
                                 verticalalignment='center',
                                 fontsize=ticklabels)
    for j in range(0, num, step):
        axs[-1, j].set_xlabel(f'{logess[j]:2.0f}',
                              fontsize=ticklabels)

    axss.append(axs)

for g, given in enumerate(givens):

    low = my.a2eq(given, AA, RR)

    high = highs[g]
    axs = axss[g]

    T = my.fitness(high, low, given, AA, RR)
    R = my.fitness(high, high, given, AA, RR)
    P = my.fitness(low, low, given, AA, RR)
    S = my.fitness(low, high, given, AA, RR)
    Z = my.gamecolors(T, R, P, S)
    greys = np.full([*Z.shape], [0.8, 0.8, 0.8, 1.0])
    mask = (Z == my.colormap['white'])
    Z[mask] = greys[mask]

    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            y = [T[i, j], R[i, j], P[i, j], S[i, j]]
            for line in axs[i, j].get_lines():
                line.remove()
            axs[i, j].plot(xaxis,
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
