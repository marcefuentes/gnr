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
alpha = 0.7
loges = 5.0
ext = 201

fslarge = 32 # Label font size
fssmall = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

rho = 1.0 - 1.0/pow(2, loges)
RRR, AAA = np.meshgrid(np.repeat(rho, ext),
                        np.repeat(alpha, ext))
xmin = 0.0
xmax = mymodule.a2max
ymin = 0.0
ymax = mymodule.a2max
x = np.linspace(xmin, xmax, num=ext)
y = np.flip(x)
X, Y = np.meshgrid(x, y)
T = mymodule.fitness(Y, X, given, AAA, RRR)
R = mymodule.fitness(Y, Y, given, AAA, RRR)
P = mymodule.fitness(X, X, given, AAA, RRR)
S = mymodule.fitness(X, Y, given, AAA, RRR)

xlabel = 'Effort to get $\it{B}$'
ylabel = 'Effort to get $\it{B}$'
letter = ord('a')
letterposition = ext*1.035
xticks = [0, ext/2, ext]
yticks = [0, ext/2, ext]
xticklabels = [f'{xmin:3.1f}',
                f'{(xmin + xmax)/2.0:3.1f}',
                f'{xmax:3.1f}']
yticklabels = [f'{ymin:3.1f}',
                f'{(ymin + ymax)/2.0:3.1f}',
                f'{ymax:3.1f}']
extent = 0, ext, 0, ext

fig, axs = plt.subplots(nrows=1,
                        ncols=len(titles),
                        figsize=(6*len(titles), 6))
fig.supxlabel(xlabel,
                x=0.513,
                y=0.0,
                fontsize=fslarge*1.2)
fig.supylabel(ylabel,
                x=0.03,
                y=0.493,
                fontsize=fslarge*1.2)

for ax, title in zip(axs, titles):
    ax.set(xticks=xticks, yticks=yticks)
    ax.set_xticklabels(xticklabels, fontsize=fssmall)
    ax.text(0,
            letterposition,
            chr(letter),
            fontsize=fslarge*0.8,
            weight='bold')
    letter += 1
    ax.set_title(title, pad=30.0, fontsize=fslarge)
    if ax.get_subplotspec().is_first_col():
        ax.set_yticklabels(yticklabels, fontsize=fssmall)
    else:
        ax.set_yticklabels([])

maskxy = (X >= Y)
N = np.full([ext, ext, 4], [1.0, 1.0, 1.0, 0.0])
masknodilemma = (mymodule.harmony(T, R, P, S) | (mymodule.deadlock(T, R, P, S) & (2.0*P > T + S)))
N[masknodilemma] = [1.0, 1.0, 1.0, 1.0]
G = np.full([ext, ext, 4], [1.0, 1.0, 1.0, 0.0])
G[maskxy] = [0.9, 0.9, 0.9, 1.0]

Z = np.full([ext, ext, 4], mymodule.colormap['white'])
mymodule.gamecolors(T, R, P, S, Z)
Z[maskxy] = [0.9, 0.9, 0.9, 1.0]
axs[0].imshow(Z, extent=extent)

Z = np.full([ext, ext], -3.0)
Z = R - P
axs[1].imshow(Z, extent=extent, vmin=-1, vmax=1)
axs[1].imshow(N, extent=extent)
axs[1].imshow(G, extent=extent)

Z = np.full([ext, ext], -3.0)
mask = mymodule.dilemma(T, R, P, S)
Z = T + S - 2.0*R
axs[2].imshow(Z, extent=extent, vmin=-1, vmax=1)
axs[2].imshow(N, extent=extent)
axs[2].imshow(G, extent=extent)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
