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
ext = 512

fslarge = 32 # Label font size
fssmall = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=num)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
xmin = 0.0
xmax = mymodule.a2max
ymin = 0.0
ymax = mymodule.a2max
X, Y = np.meshgrid(np.linspace(xmin, xmax, num=ext),
                    np.linspace(ymax, ymin, num=ext))

cmap = plt.cm.viridis
cmap.set_bad(color='white')
step = int(num/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
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
for outer, title in zip(outergrid, titles):
    grid = outer.subgridspec(nrows=num,
                                ncols=num,
                                wspace=0,
                                hspace=0)
    axs = grid.subplots()
    axs[0, int(num/2)].set_title(f'{title}',
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
    if title == 'Games':
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
    AAA = np.full([ext, ext], alpha)
    for j, rho in enumerate(rhos):
        RRR = np.full([ext, ext], rho)
        T = mymodule.fitness(Y, X, given, AAA, RRR)
        R = mymodule.fitness(Y, Y, given, AAA, RRR)
        P = mymodule.fitness(X, X, given, AAA, RRR)
        S = mymodule.fitness(X, Y, given, AAA, RRR)
        Z = np.full([ext, ext, 4], mymodule.colormap['white'])
        mymodule.gamecolors(T, R, P, S, Z)
        mask = (X > Y)
        Z[mask] = [0.9, 0.9, 0.9, 1.0]
        axss[0][i][j].imshow(Z, extent=extent)
        #axs[i, j].plot(0.01*a2eq[i, j]*ext,
        #            (0.99*a2social[i, j] + 0.01*mymodule.a2max)*ext,
        #            marker='o',
        #            color='orange',
        #            markersize=4)
        #axs[i, j].plot(0.99*a2eq[i, j]*ext,
        #            (0.01*a2social[i, j] + 0.99*mymodule.a2max)*ext,
        #            marker='o',
        #            color='red',
        #            markersize=4)

        Z = np.zeros([ext, ext])
        mask = mymodule.dilemma(T, R, P, S)
        Z[mask] = R[mask] - P[mask] + 0.000001
        Z = np.ma.masked_where(Z == 0.0, Z)
        axss[1][i][j].imshow(Z, extent=extent, cmap=cmap, vmin=-1, vmax=1)

        Z = np.zeros([ext, ext])
        mask = mymodule.dilemma(T, R, P, S)
        Z[mask] = 1.0 - (2.0*R[mask] - T[mask] - S[mask])
        Z = np.ma.masked_where(Z == 0.0, Z)
        axss[2][i][j].imshow(Z, extent=extent, cmap=cmap, vmin=-1, vmax=1)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
