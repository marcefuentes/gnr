#! /usr/bin/env python

import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter ()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

givens = [0.95]
#givens = np.linspace(0.0, 1.0, num=21)

num = 21    # Number of subplot rows and columns

fslarge = 32 # Label font size
fssmall = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

if givens[-1] > 0.9999999:
    givens[-1] = 0.9999999
alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=num)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
a2social = mymodule.a2eq(0.0, AA, RR)
highs = [np.full([num, num], 0.5*a2social + 0.5*mymodule.a2max), 
        np.full([num, num], 0.01*a2social + 0.99*mymodule.a2max)] 

xlim=[0, 5]
ylim=[0.0, 2.0]
step = int(num/2)
xaxis = [1, 2, 3, 4]
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'

fig = plt.figure(figsize=(16, 8))
fig.supxlabel(xlabel,
                x=0.525,
                y=0.03,
                fontsize=fslarge)
fig.supylabel(ylabel,
                x=0.05,
                y=0.52,
                fontsize=fslarge)

outergrid = fig.add_gridspec(nrows=1,
                                ncols=len(highs),
                                left=0.15,
                                right=0.9,
                                top=0.86,
                                bottom=0.176)
axss = []
for l, outer in enumerate(outergrid):
    grid = outer.subgridspec(nrows=num,
                                ncols=num,
                                wspace=0,
                                hspace=0)
    axs = grid.subplots()
    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            axs[i, j].set(xticks=[],
                            yticks=[],
                            xlim=xlim,
                            ylim=ylim)
    if l == 0:
        for i in range(0, num, step):
            axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                                    rotation='horizontal',
                                    horizontalalignment='right',
                                    verticalalignment='center',
                                    fontsize=fssmall)
    for j in range(0, num, step):
        axs[-1, j].set_xlabel(f'{logess[j]:2.0f}', fontsize=fssmall)

    axss.append(axs)

frames = []
for given in givens:

    a2eq = mymodule.a2eq(given, AA, RR)
    lows = [np.full([num, num], 0.5*a2eq),
            np.full([num, num], 0.99*a2eq)]

    for axs, low, high in zip(axss, lows, highs): 
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
                                markersize=3)

    text = fig.text(0.90,
                    0.043,
                    'Given: ' + f'{given:4.2f}',
                    fontsize=fslarge,
                    color='grey',
                    ha='right')
    plt.savefig('temp.png', transparent=False)
    text.remove()
    frames.append(iio.imread('temp.png'))
    os.remove('temp.png')

plt.close()

iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
