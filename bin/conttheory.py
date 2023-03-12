#! /usr/bin/env python

import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule as my
import numpy as np
import os
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

titles = ['Games',
          '$\it{R}$ - $\it{P}$',
          '$\it{T}$ + $\it{S}$ - 2$\it{R}$']
givens = np.linspace(0.95, 1.0, num=1)
#givens = np.linspace(0.0, 1.0, num=21)

num = 21    # Number of subplot rows & columns
ext = 256

plotsize = 20

alphas = np.linspace(my.alphamax, my.alphamin, num=num)
logess = np.linspace(my.logesmin, my.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)

step = int(num/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
letterposition = ext*1.035
extent = 0, ext, 7.5, ext
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
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(0.1)
    if g == 0:
        for i in range(0, num, step):
            axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                                 rotation='horizontal',
                                 horizontalalignment='right',
                                 verticalalignment='center',
                                 fontsize=ticklabels)
    for j in range(0, num, step):
        axs[-1, j].set_xlabel(f'{logess[j]:2.0f}',
                              x=0.45,
                              fontsize=ticklabels)

    axss.append(axs)

frames = []
for given in givens:

    for i, alpha in enumerate(alphas):
        AA = np.full([ext, ext], alpha)
        for j, rho in enumerate(rhos):

            xmin = 0.0
            xmax = my.a2eq(given, alpha, rho)
            ymin = my.a2eq(0.0, alpha, rho)
            ymax = my.a2max
            x = np.linspace(xmin, xmax, num=ext)
            y = np.linspace(ymax, ymin, num=ext)
            X, Y = np.meshgrid(x, y)
            G = np.full([ext, ext, 4], my.colormap['transparent'])
            G[X >= Y] = [0.9, 0.9, 0.9, 1.0]
            RR = np.full([ext, ext], rho)
            T = my.fitness(Y, X, given, AA, RR)
            R = my.fitness(Y, Y, given, AA, RR)
            P = my.fitness(X, X, given, AA, RR)
            S = my.fitness(X, Y, given, AA, RR)

            Z = my.gamecolors(T, R, P, S)
            Z[X >= Y] = [0.9, 0.9, 0.9, 1.0]
            axss[0][i][j].imshow(Z, extent=extent)

            Z = R - P
            axss[1][i][j].imshow(Z, extent=extent, vmin=-1, vmax=1)
            axss[1][i][j].imshow(G, extent=extent)

            Z = T + S - 2.0*R
            m = R < P
            Z[m] = T[m] + S[m] - 2.0*P[m]
            axss[2][i][j].imshow(Z, extent=extent, vmin=-1, vmax=1)
            axss[2][i][j].imshow(G, extent=extent)

    text = fig.text(0.85,
                    0.02,
                    'Given: ' + f'{given:4.2f}',
                    fontsize=biglabels,
                    color='white',
                    ha='right')
    plt.savefig('temp.png', transparent=False)
    text.remove()
    frames.append(iio.imread('temp.png'))
    os.remove('temp.png')

plt.close()

iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
