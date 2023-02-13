#! /usr/bin/env python

import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import time

start_time = time.perf_counter ()

givenmin = 0.95
givenmax = 0.95

num = 21    # Number of subplot rows and columns
ngiven = 21
filename = 'trps'

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

every = int(num/2)

if givenmin != givenmax:
    movie = True
    givens = np.linspace(givenmin, givenmax, num=ngiven)
    frames = []
else:
    movie = False 
    givens = np.array([givenmin])

nc = num
nr = num
alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=nr)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=nc)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)

xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'

zeros = np.zeros([nr, nc])
low = np.full([nr, nc], mymodule.a2low)
high = np.full([nr, nc], mymodule.a2high)
R = mymodule.fitness(high, high, zeros, AA, RR)
P = mymodule.fitness(low, low, zeros, AA, RR)
xaxis = [1, 2, 3, 4]

fig = plt.figure(figsize=(8, 8))
fig.supxlabel(xlabel, x=0.56, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.05, y=0.52, fontsize=fslabel)
grid = fig.add_gridspec(nrows=nr, ncols=nc, left=0.22, right=0.9, top=0.86, bottom=0.176, wspace=0, hspace=0)
axs = grid.subplots()

for row in axs:
    for ax in row:
        ax.set(xticks=[], yticks=[], xlim=(0, 5), ylim=(0.0, 2.0))
for ax, loges in zip(axs[-1, ::every], logess[::every]):
    ax.set_xlabel(round(loges), fontsize=fstick)
for ax, alpha in zip(axs[::every, 0], alphas[::every]):
    ax.set_ylabel(f'{alpha:1.1f}', rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

for given in givens:

    T = mymodule.fitness(high, low, given, AA, RR)
    S = mymodule.fitness(low, high, given, AA, RR)
    Z = np.full([nr, nc, 4], mymodule.colormap['default'])
    mymodule.gametypes(T, R, P, S, Z)

    for row, rowT, rowR, rowP, rowS, rowZ in zip(axs, T, R, P, S, Z):
        for ax, tt, rr, pp, ss, zz in zip(row, rowT, rowR, rowP, rowS, rowZ):
            y = [tt, rr, pp, ss]
            for line in ax.get_lines():
                line.remove()
            ax.plot(xaxis, y, c=zz, linewidth=3, marker='o', markerfacecolor='white', markersize=3)

    text = fig.text(0.90,
                    0.035,
                    f'given\n{given:4.2f}',
                    fontsize=fstick+4,
                    color='grey',
                    ha='right')

    if movie:
        plt.savefig('temp.png', transparent=False)
        text.remove()
        frames.append(iio.imread('temp.png'))
        os.remove('temp.png')
    else:
        plt.savefig(filename + '.png', transparent=False)

plt.close()

if movie:
    iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
