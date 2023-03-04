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

titles = ['Games',
          '$\it{R}$ - $\it{P}$',
          '$\it{T}$ + $\it{S}$ - 2$\it{R}$']
#givens = np.linspace(0.0, 1.0, num=21)
givens = np.linspace(0.95, 1.0, num=1)
alpha = 0.46
loges = 2.5
ext = 1024

plotsize = 6

rho = 1.0 - 1.0/pow(2, loges)
RRR, AAA = np.meshgrid(np.repeat(rho, ext),
                        np.repeat(alpha, ext))
x = np.linspace(0.0, mymodule.a2max, num=ext)
y = np.flip(x)
X, Y = np.meshgrid(x, y)
G = np.full([ext, ext, 4], mymodule.colormap['transparent'])
maskxy = (X >= Y)
G[maskxy] = [0.9, 0.9, 0.9, 1.0]

xlabel = 'Effort to get $\it{B}$'
ylabel = 'Effort to get $\it{B}$'
letter = ord('a')
letterposition = 1.035
xticks = [-0.5, ext/2-0.5, ext-0.5]
yticks = [-0.5, ext/2-0.5, ext-0.5]
xticklabels = [f'{0.0:3.1f}',
               f'{mymodule.a2max/2.0:3.1f}',
               f'{mymodule.a2max:3.1f}']
yticklabels = np.flip(xticklabels)
width = plotsize*len(titles)
height = plotsize
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig, axs = plt.subplots(nrows=1,
                        ncols=len(titles),
                        figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.513,
              y=0.0,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.04,
              y=0.493,
              fontsize=biglabels)

for ax in fig.get_axes():
    ax.set(xticks=xticks, yticks=yticks)
    ax.set_xticklabels(xticklabels, fontsize=ticklabels)
    ax.set_yticklabels([])
    ax.text(0,
            letterposition,
            chr(letter),
            transform=ax.transAxes,
            fontsize=plotsize*5,
            weight='bold')
    letter += 1
for j, title in enumerate(titles):
    axs[j].set_title(title, pad=plotsize*5, fontsize=plotsize*5)
axs[0].set_yticklabels(yticklabels, fontsize=ticklabels)

frames = []
for given in givens:

    T = mymodule.fitness(Y, X, given, AAA, RRR)
    R = mymodule.fitness(Y, Y, given, AAA, RRR)
    P = mymodule.fitness(X, X, given, AAA, RRR)
    S = mymodule.fitness(X, Y, given, AAA, RRR)

    Z = mymodule.gamecolors(T, R, P, S)
    Z[maskxy] = [0.9, 0.9, 0.9, 1.0]
    axs[0].imshow(Z)

    N = mymodule.nodilemmacolors(T, R, P, S)

    Z = R - P
    axs[1].imshow(Z, vmin=-1, vmax=1)
    axs[1].imshow(N)
    axs[1].imshow(G)

    Z = T + S - 2.0*R
    axs[2].imshow(Z, vmin=-1, vmax=1)
    axs[2].imshow(N)
    axs[2].imshow(G)

    text = fig.text(0.90,
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
