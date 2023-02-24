#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter ()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

titles = ['Games',
                '$\it{R}$ - $\it{P}$',
                '$\it{T}$ + $\it{S}$ - 2$\it{R}$']
given = 0.95
num = 5    # Number of subplot rows & columns
ext = 256

fslarge = 32 # Label font size
fssmall = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

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
G = np.full([ext, ext, 4], [1.0, 1.0, 1.0, 0.0])
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
xticks = [0, ext/2, ext]
yticks = [0, ext/2, ext]
xticklabels = [f'{xmin:2.0f}',
                f'{(xmin + xmax)/2.0:2.0f}',
                f'{xmax:2.0f}']
yticklabels = [f'{ymin:3.1f}',
                f'{(ymin + ymax)/2.0:3.1f}',
                f'{ymax:3.1f}']
extent = 0, ext, 0, ext

fig = plt.figure(figsize=(6*len(titles), 6))
fig.supxlabel(xlabel,
                x=0.525,
                y=0.0,
                fontsize=fslarge)
fig.supylabel(ylabel,
                x=0.08,
                y=0.52,
                fontsize=fslarge)

outergrid = fig.add_gridspec(nrows=1,
                                ncols=len(titles),
                                left=0.15,
                                right=0.9,
                                top=0.86,
                                bottom=0.176)

axss = []
for g, title in enumerate(titles):
    grid = outergrid[g].subgridspec(nrows=num,
                                    ncols=num,
                                    wspace=0,
                                    hspace=0)
    axs = grid.subplots()
    axs[0, int(num/2)].set_title(title,
                                    pad=30.0,
                                    fontsize=fslarge*0.8)
    axs[0, 0].set_title(chr(letter),
                        fontsize=fslarge*0.8,
                        weight='bold',
                        loc='left')
    letter += 1

    for ax in fig.get_axes():
        ax.set(xticks=[],
                yticks=[])
    if g == 0:
        for i in range(0, num, step):
            axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                                    rotation='horizontal',
                                    horizontalalignment='right',
                                    verticalalignment='center',
                                    fontsize=fssmall)
    for j in range(0, num, step):
        axs[-1, j].set_xlabel(f'{logess[j]:2.0f}', fontsize=fssmall)

    axss.append(axs)

for i, alpha in enumerate(alphas):
    AA = np.full([ext, ext], alpha)
    for j, rho in enumerate(rhos):
        RR = np.full([ext, ext], rho)
        T = mymodule.fitness(Y, X, given, AA, RR)
        R = mymodule.fitness(Y, Y, given, AA, RR)
        P = mymodule.fitness(X, X, given, AA, RR)
        S = mymodule.fitness(X, Y, given, AA, RR)

        N = np.full([ext, ext, 4], [1.0, 1.0, 1.0, 0.0])
        masknodilemma = (mymodule.harmony(T, R, P, S) | (mymodule.deadlock(T, R, P, S) & (2.0*P > T + S)))
        N[masknodilemma] = [1.0, 1.0, 1.0, 1.0]

        Z = np.full([ext, ext, 4], mymodule.colormap['white'])
        mymodule.gamecolors(T, R, P, S, Z)
        Z[maskxy] = [0.9, 0.9, 0.9, 1.0]
        axss[0][i][j].imshow(Z, extent=extent)

        Z = np.full([ext, ext], -3.0)
        Z = R - P
        axss[1][i][j].imshow(Z, extent=extent, vmin=-1, vmax=1)
        axss[1][i][j].imshow(N, extent=extent)
        axss[1][i][j].imshow(G, extent=extent)

        Z = np.full([ext, ext], -3.0)
        mask = mymodule.dilemma(T, R, P, S)
        Z = T + S - 2.0*R
        axss[2][i][j].imshow(Z, extent=extent, vmin=-1, vmax=1)
        axss[2][i][j].imshow(N, extent=extent)
        axss[2][i][j].imshow(G, extent=extent)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
