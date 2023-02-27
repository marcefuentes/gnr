#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

titles = ['Games',
          '$\it{R}$ - $\it{P}$',
          '$\it{T}$ + $\it{S}$ - 2$\it{R}$']
given = 0.95
num = 5    # Number of subplot rows & columns
ext = 256

plotsize = 6

alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=num)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
xmin = 0.0
xmax = mymodule.a2max
ymin = 0.0
ymax = mymodule.a2max
x = np.linspace(xmin, xmax, num=ext)
y = np.flip(x)
X, Y = np.meshgrid(x, y)
G = np.full([ext, ext, 4], mymodule.colormap['transparent'])
maskxy = (X >= Y)
G[maskxy] = [0.9, 0.9, 0.9, 1.0]

step = int(num/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
letterposition = ext*1.035
xmin = logess[0]
xmax = logess[-1]
ymin = alphas[-1]
ymax = alphas[0]
extent = 0, ext, 0+7.5, ext
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
    if g == 0:
        for i in range(0, num, step):
            axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                                 rotation='horizontal',
                                 horizontalalignment='right',
                                 verticalalignment='center',
                                 fontsize=ticklabels)
    for j in range(0, num, step):
        axs[-1, j].set_xlabel(f'{logess[j]:2.0f}', fontsize=ticklabels)

    axss.append(axs)

for i, alpha in enumerate(alphas):
    AA = np.full([ext, ext], alpha)
    for j, rho in enumerate(rhos):
        RR = np.full([ext, ext], rho)
        T = mymodule.fitness(Y, X, given, AA, RR)
        R = mymodule.fitness(Y, Y, given, AA, RR)
        P = mymodule.fitness(X, X, given, AA, RR)
        S = mymodule.fitness(X, Y, given, AA, RR)

        Z = mymodule.gamecolors(T, R, P, S)
        Z[maskxy] = [0.9, 0.9, 0.9, 1.0]
        axss[0][i][j].imshow(Z, extent=extent)

        N = mymodule.nodilemmacolors(T, R, P, S)

        Z = R - P
        axss[1][i][j].imshow(Z, extent=extent, vmin=-1, vmax=1)
        axss[1][i][j].imshow(N, extent=extent)
        axss[1][i][j].imshow(G, extent=extent)

        Z = T + S - 2.0*R
        axss[2][i][j].imshow(Z, extent=extent, vmin=-1, vmax=1)
        axss[2][i][j].imshow(N, extent=extent)
        axss[2][i][j].imshow(G, extent=extent)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
