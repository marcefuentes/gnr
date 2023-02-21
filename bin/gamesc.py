#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter ()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

givens = [0.5, 0.95]
num = 5    # Number of subplot rows & columns
numa2 = 512

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
xmin = 0.0
xmax = mymodule.a2max
ymin = 0.0
ymax = mymodule.a2max

X, Y = np.meshgrid(np.linspace(xmin, xmax, num=numa2),
                    np.linspace(ymax, ymin, num=numa2))

every = int(num/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
extent2 = 0, numa2, 0, numa2

fig = plt.figure(figsize=(6*len(givens), 6))
fig.supxlabel(xlabel,
                x=0.525,
                y=0.0,
                fontsize=fslarge)
fig.supylabel(ylabel,
                x=0.05,
                y=0.52,
                fontsize=fslarge)

outergrid = fig.add_gridspec(nrows=1,
                                ncols=len(givens),
                                left=0.15,
                                right=0.9,
                                top=0.86,
                                bottom=0.176)

for outer, given in zip(outergrid, givens):
    grid = outer.subgridspec(nrows=num,
                                ncols=num,
                                wspace=0,
                                hspace=0)
    axs = grid.subplots()
    axs[0, int(num/2)].set_title(f'{given*100:2.0f}%',
                                    pad=30.0,
                                    fontsize=fslarge*0.8)
    axs[0, 0].set_title(chr(letter),
                        fontsize=fslarge*0.8,
                        weight='bold',
                        loc='left')
    letter += 1

    for i in range(num):
        for j in range(num):
            axs[i, j].set(xticks=[],
                    yticks=[])
    for ax, loges in zip(axs[-1, ::every], logess[::every]):
        ax.set_xlabel(f'{loges:2.0f}', fontsize=fssmall)
    if given == 0.5:
        for ax, alpha in zip(axs[::every, 0], alphas[::every]):
            ax.set_ylabel(f'{alpha:3.1f}',
                            rotation='horizontal',
                            horizontalalignment='right',
                            verticalalignment='center',
                            fontsize=fssmall)

    MRT = MRT0*(1.0 - given)
    Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)

    for i in range(num): 
        for j in range(num):
            AAA = np.full([numa2, numa2], alphas[i])
            RRR = np.full([numa2, numa2], rhos[j])
            T = mymodule.fitness(Y, X, given, AAA, RRR)
            R = mymodule.fitness(Y, Y, given, AAA, RRR)
            P = mymodule.fitness(X, X, given, AAA, RRR)
            S = mymodule.fitness(X, Y, given, AAA, RRR)
            Z = np.full([numa2, numa2, 4], mymodule.colormap['white'])
            mymodule.gamecolors(T, R, P, S, Z)
            mask = (X > Y)
            Z[mask] = [0.9, 0.9, 0.9, 1.0]
            axs[i, j].imshow(Z, extent=extent2)
            axs[i, j].plot(0.01*a2eq[i, j]*numa2,
                        (0.99*a2social[i, j] + 0.01*mymodule.a2max)*numa2,
                        marker='o',
                        color='orange',
                        markersize=4)
            axs[i, j].plot(0.99*a2eq[i, j]*numa2,
                        (0.01*a2social[i, j] + 0.99*mymodule.a2max)*numa2,
                        marker='o',
                        color='red',
                        markersize=4)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
