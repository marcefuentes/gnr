#! /usr/bin/env python

import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import time

start_time = time.perf_counter ()

traitlabels = ['Effort to get $\it{B}$', 'Fitness']
alphamin = 0.1
alphamax = 0.9
logesmin = -5.0
logesmax = 5.0
givenmin = 0.95
givenmax = 0.95

num = 21    # Number of subplot rows and columns
ngiven = 21
filename = 'gamesd'

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

if givenmin != givenmax:
    movie = True
    givens = np.linspace(givenmin, givenmax, num=ngiven)
    frames = []
else:
    movie = False 
    givens = np.array([givenmin])

nc = num
nr = num
alphas = np.linspace(alphamax, alphamin, num=nr)
logess = np.linspace(logesmin, logesmax, num=nc)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)

xmin = logesmin
xmax = logesmax
xlabel = 'Substitutability of $\it{B}$'
ymin = alphamin
ymax = alphamax
ylabel = 'Value of $\it{B}$'

traitvmaxs = [mymodule.a2max,
                mymodule.fitness(np.array([mymodule.a2max]),
                np.array([mymodule.a2max]),
                np.array([0.0]),
                np.array([0.9]),
                np.array([5.0]))]
xticklabels = [round(xmin), round((xmin + xmax)/2), round(xmax)]
yticklabels = [round(ymin, 1), round((ymin + ymax)/2, 1), round(ymax, 1)]
extent = 0, nc, 0, nr

zeros = np.zeros([nr, nc])
low = np.full([nr, nc], mymodule.a2low)
high = np.full([nr, nc], mymodule.a2high)
R = mymodule.fitness(high, high, zeros, AA, RR)
P = mymodule.fitness(low, low, zeros, AA, RR)
a2social = np.copy(low)
mask = (R > P)
a2social[mask] = high[mask]
wsocial = mymodule.fitness(a2social, a2social, zeros, AA, RR)

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 18))
fig.supxlabel(xlabel, x=0.515, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.04, y=0.510, fontsize=fslabel)
cmap = plt.cm.viridis

letter = ord('a')
for axrow in axs:
    for ax in axrow:
        ax.text(0,
                nr*1.035,
                chr(letter),
                fontsize=fslabel,
                weight='bold')
        letter += 1
        ax.set(xticks=[0, nc/2, nc],
                yticks=[0, nr/2, nr],
                xticklabels=[],
                yticklabels=[])
        if ax.get_subplotspec().is_first_row():
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

for given in givens:

    T = mymodule.fitness(high, low, given, AA, RR)
    S = mymodule.fitness(low, high, given, AA, RR)

    Z = np.full([nr, nc, 4], mymodule.colormap['default'])
    TS = np.copy(zeros)
    mymodule.gametypes(T, R, P, S, Z, TS)
    ax = axs[0, 0]
    ax.imshow(Z, extent=extent)
    ax.set_title('Game types', pad=50.0, fontsize=fslabel)
    ax = axs[0, 1]
    ax.imshow(TS, extent=extent, cmap=cmap, vmin=0, vmax=0.7)
    ax.set_title('Prisoner\'s dilemma\n$\it{T}$ + $\it{S}$ - 2$\it{R}$',
                    pad=50.0,
                    fontsize=fslabel)

    a2eq = np.copy(zeros)
    weq = np.copy(zeros)
    mymodule.equilibrium(T, R, P, S, low, high, a2eq, weq)
    Mss = [[a2social, wsocial], [a2eq, weq]]

    for row, Ms in zip(axs[1:], Mss):
        for ax, M, traitvmax in zip(row, Ms, traitvmaxs):
            ax.imshow(M, extent=extent, cmap=cmap, vmin=0, vmax=traitvmax)

    if movie:
        text = fig.text(0.90,
                        0.93,
                        f'given\n{given:4.2f}',
                        fontsize=fstick+4,
                        color='grey',
                        ha='right')
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
