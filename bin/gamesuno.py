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
givens = np.linspace(0.0, 1.0, num=21)

given = 0.95
alpha = 0.7
loges = 5.0
ext = 201

fslarge = 32 # Label font size
fssmall = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

if givens[-1] > 0.9999999:
    givens[-1] = 0.9999999
rho = 1.0 - 1.0/pow(2, loges)
#MRT0 = mymodule.b*mymodule.Rq
#Q0 = mymodule.Rq*pow(MRT0*alpha/(1.0 - alpha), 1.0/(rho - 1.0))
#a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)
#MRT = MRT0*(1.0 - given)
#Q = mymodule.Rq*pow(MRT*alpha/(1.0 - alpha), 1.0/(rho - 1.0))
#a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
xmin = 0.0
xmax = mymodule.a2max
ymin = 0.0
ymax = mymodule.a2max
x = np.linspace(xmin, xmax, num=ext)
y = np.flip(x)
X, Y = np.meshgrid(x, y)
RRR, AAA = np.meshgrid(np.repeat(rho, ext),
                        np.repeat(alpha, ext))
T = mymodule.fitness(Y, X, given, AAA, RRR)
R = mymodule.fitness(Y, Y, given, AAA, RRR)
P = mymodule.fitness(X, X, given, AAA, RRR)
S = mymodule.fitness(X, Y, given, AAA, RRR)

cmap = plt.cm.viridis
cmap.set_bad(color='white')
xlabel = 'Effort to get $\it{B}$'
ylabel = 'Effort to get $\it{B}$'
xticks = [0, ext/2, ext]
yticks = [0, ext/2, ext]
xticklabels = [f'{xmin:3.1f}',
                f'{(xmin + xmax)/2.0:3.1f}',
                f'{xmax:3.1f}']
yticklabels = [f'{ymin:3.1f}',
                f'{(ymin + ymax)/2.0:3.1f}',
                f'{ymax:3.1f}']
letter = ord('a')
letterposition = ext*1.035
extent = 0, ext, 0, ext

fig, axs = plt.subplots(nrows=1,
                        ncols=len(titles),
                        figsize=(6*len(titles), 6))
fig.supxlabel(xlabel,
                x=0.513,
                y=0.01,
                fontsize=fslarge*1.2)
fig.supylabel(ylabel,
                x=0.03,
                y=0.493,
                fontsize=fslarge*1.2)

for i, title in enumerate(titles):
    ax = axs[i]
    ax.text(0, 
            letterposition,
            chr(letter),
            fontsize=fslarge*0.8,
            weight='bold')
    letter += 1
    ax.set(xticks=xticks, yticks=yticks)
    if ax.get_subplotspec().is_first_col():
        ax.set_yticklabels(yticklabels, fontsize=fssmall)
    else:
        ax.set_yticklabels([])
    ax.set_title(title, pad=40.0, fontsize=fslarge*0.9)
    ax.set_xticklabels(xticklabels, fontsize=fssmall)

maskxy = (X >= Y) 
G = np.full([ext, ext, 4], [0.0, 0.0, 0.0, 0.0])
G[maskxy] = [0.9, 0.9, 0.9, 1.0]

Z = np.full([ext, ext, 4], mymodule.colormap['white'])
mymodule.gamecolors(T, R, P, S, Z)
Z[maskxy] = [0.7, 0.7, 0.7, 1.0]
axs[0].imshow(Z, extent=extent)

Z = np.full([ext, ext], -3.0)
mask = mymodule.dilemma(T, R, P, S)
Z[mask] = R[mask] - P[mask]
Z = np.ma.masked_where(Z == -3.0, Z)
axs[1].imshow(Z, extent=extent, cmap=cmap, vmin=-1, vmax=1)
axs[1].imshow(G, extent=extent)

Z = np.full([ext, ext], -3.0)
mask = mymodule.dilemma(T, R, P, S)
Z[mask] = 1.0 - (2.0*R[mask] - T[mask] - S[mask])
Z = np.ma.masked_where(Z == -3.0, Z)
axs[2].imshow(Z, extent=extent, cmap=cmap, vmin=-1, vmax=1)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
