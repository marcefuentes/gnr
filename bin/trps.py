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
MRT0 = mymodule.b*mymodule.Rq
Q0 = mymodule.Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)

xlim=[0, 5]
ylim=[0.0, 2.0]
every = int(num/2)
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
                                ncols=2,
                                left=0.15,
                                right=0.9,
                                top=0.86,
                                bottom=0.176)
axss = []
for outer in outergrid:
    grid = outer.subgridspec(nrows=num,
                                ncols=num,
                                wspace=0,
                                hspace=0)
    axss.append(grid.subplots())

for axs in axss:
    for axrow in axs:
        for ax in axrow:
            ax.set(xticks=[],
                    yticks=[],
                    xlim=xlim,
                    ylim=ylim)
    for ax, loges in zip(axs[-1, ::every], logess[::every]):
        ax.set_xlabel(round(loges), fontsize=fssmall)
for ax, alpha in zip(axss[0][::every, 0], alphas[::every]):
    ax.set_ylabel(f'{alpha:1.1f}',
                    rotation='horizontal',
                    horizontalalignment='right',
                    verticalalignment='center',
                    fontsize=fssmall)

frames = []
for given in givens:

    MRT = MRT0*(1.0 - given)
    Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
    a2lows = [0.01*a2eq, 0.99*a2eq]
    a2highs = [(99.0*a2social + mymodule.a2max)/100.0,
                (a2social + 99.0*mymodule.a2max)/100.0]

    for axs, a2low, a2high in zip(axss, a2lows, a2highs): 
        low = np.full([num, num], a2low)
        high = np.full([num, num], a2high)
        T = mymodule.fitness(high, low, given, AA, RR)
        R = mymodule.fitness(high, high, given, AA, RR)
        P = mymodule.fitness(low, low, given, AA, RR)
        S = mymodule.fitness(low, high, given, AA, RR)
        Z = np.full([num, num, 4], mymodule.colormap['red'])
        mymodule.gamecolors(T, R, P, S, Z)

        for i in range(num):
            for j in range(num):
                y = [T[i, j], R[i, j], P[i, j], S[i, j]]
                for line in ax.get_lines():
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
