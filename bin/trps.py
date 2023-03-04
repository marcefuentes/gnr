#! /usr/bin/env python

import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

givens = np.linspace(0.95, 1.0, num=1)
#givens = np.linspace(0.0, 1.0, num=21)
distances = np.linspace(0.8, 0.2, num=3)
titles = []
for distance in distances:
    titles.append(f'{distance*100:2.0f}%')

num = 21    # Number of subplot rows and columns
plotsize = 8

alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=num)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
highs = [] 
eq = mymodule.a2eq(0.0, AA, RR)
for distance in distances:
    highs.append((1.0 - distance)*eq + distance*mymodule.a2max)

xlim=[0, 5]
ylim=[0.0, 2.0]
step = int(num/2)
xaxis = [1, 2, 3, 4]
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
width = plotsize*len(titles)
height = plotsize
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.502,
              y=0.0,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.07,
              y=0.5,
              fontsize=biglabels)

outergrid = fig.add_gridspec(nrows=1,
                             ncols=len(titles),
                             left=0.15,
                             right=0.85,
                             top=0.8,
                             bottom=0.2)

axss = []
for g, title in enumerate(titles):
    grid = outergrid[g].subgridspec(nrows=num,
                                    ncols=num,
                                    wspace=0,
                                    hspace=0)
    axs = grid.subplots()
    axs[0, int(num/2)].set_title(title,
                                 pad=plotsize*5,
                                 fontsize=plotsize*5)
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

frames = []
for given in givens:

    eq = mymodule.a2eq(given, AA, RR)
    for g, distance in enumerate(distances):

        low = distance*eq
        high = highs[g]
        axs = axss[g]

        T = mymodule.fitness(high, low, given, AA, RR)
        R = mymodule.fitness(high, high, given, AA, RR)
        P = mymodule.fitness(low, low, given, AA, RR)
        S = mymodule.fitness(low, high, given, AA, RR)
        Z = mymodule.gamecolors(T, R, P, S)

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

    text = fig.text(0.85,
                    0.02,
                    'Given: ' + f'{given:4.2f}',
                    fontsize=biglabels,
                    color='grey',
                    ha='right')
    plt.savefig('temp.png', transparent=False)
    text.remove()
    frames.append(iio.imread('temp.png'))
    os.remove('temp.png')

plt.close()

iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
