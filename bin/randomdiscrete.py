#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule as my
import numpy as np
import os
import random
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

titles = ['Games',
          '$\it{R}$ - $\it{P}$',
          '$\it{T}$ + $\it{S}$ - 2$\it{R}$']
givens = [1.0, 0.95, 0.50]
rows = givens
ext = 21
plotsize = 4

alphas = np.linspace(my.alphamax, my.alphamin, num=ext)
logess = np.linspace(my.logesmin, my.logesmax, num=ext)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
random.seed()
one = random.random()
two = random.random()
aBlow = min(one, two)
aBhigh = max(one, two)
f = open('aB.glo', 'w')
f.write(f'aBlow,{aBlow}\n')
f.write(f'aBhigh,{aBhigh}\n')
f.close()

low = np.full([*RR.shape], aBlow)
high = np.full([*RR.shape], aBhigh)

xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
letterposition = 1.035
xticks = [-0.5, ext/2-0.5, ext-0.5]
yticks = [-0.5, ext/2-0.5, ext-0.5]
xmin = logess[0]
xmax = logess[-1]
ymin = alphas[-1]
ymax = alphas[0]
xticklabels = [f'{xmin:2.0f}',
               f'{(xmin + xmax)/2.0:2.0f}',
               f'{xmax:2.0f}']
yticklabels = [f'{ymax:3.1f}',
               f'{(ymin + ymax)/2.0:3.1f}',
               f'{ymin:3.1f}']
width = plotsize*len(titles)
height = plotsize*len(rows)
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig, axs = plt.subplots(nrows=len(rows),
                        ncols=len(titles),
                        figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.513,
              y=0.01,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.03,
              y=0.493,
              fontsize=biglabels)
fig.text(0.90,
         0.02,
         f'{aBlow:4.2f}, {aBhigh:4.2f}',
         fontsize=biglabels,
         color='grey',
         ha='right')

for ax in fig.get_axes():
    ax.set(xticks=xticks, yticks=yticks)
    ax.set(xticklabels=[], yticklabels=[])
    ax.text(0,
            letterposition,
            chr(letter),
            transform=ax.transAxes,
            fontsize=plotsize*5,
            weight='bold')
    letter += 1
for i, row in enumerate(rows):
    axs[i, 0].set_yticklabels(yticklabels, fontsize=ticklabels)
for j, title in enumerate(titles):
    axs[0, j].set_title(title, pad=plotsize*10, fontsize=plotsize*5)
    axs[-1, j].set_xticklabels(xticklabels, fontsize=ticklabels)

for i, given in enumerate(givens):

    T = my.fitness(high, low, given, AA, RR)
    R = my.fitness(high, high, given, AA, RR)
    P = my.fitness(low, low, given, AA, RR)
    S = my.fitness(low, high, given, AA, RR)

    Z = my.gamecolors(T, R, P, S)
    axs[i, 0].imshow(Z)

    N = my.nodilemmacolors(T, R, P, S)

    Z = R - P
    axs[i, 1].imshow(Z, vmin=-1, vmax=1)
    axs[i, 1].imshow(N)

    Z = T + S - 2.0*R
    m = R < P
    Z[m] = T[m] + S[m] - 2.0*P[m]
    axs[i, 2].imshow(Z, vmin=-1, vmax=1)
    axs[i, 2].imshow(N)

plt.savefig(filename + '.png', transparent=False)
plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
