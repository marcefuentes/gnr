#! /usr/bin/env python

from math import log
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
givenmin = 0.0
givenmax = 1.0
a2low = 0.1
a2high = 0.6

num = 21    # Number of subplot rows and columns
nloges = 21
filename = 'gamesdrho2'

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

if logesmin != logesmax:
    movie = True
    logess = np.linspace(logesmin, logesmax, num=nloges)
    frames = []
else:
    movie = False 
    logess = np.array([logesmin])

nc = num
nr = num
alphas = np.linspace(alphamin, alphamax, num=nr)
givens = np.linspace(givenmax, givenmin, num=nc)
rhos = 1.0 - 1.0/pow(2, logess)
AA, GG = np.meshgrid(alphas, givens)

xmin = alphamin
xmax = alphamax
xlabel = 'Value of $\it{B}$'
ymin = givenmin
ymax = givenmax
ylabel = 'Partner\'s share of $\it{B}$'

traitvmaxs = [mymodule.a2max, mymodule.fitness(np.array([mymodule.a2max]), np.array([mymodule.a2max]), np.array([0.0]), np.array([0.9]), np.array([5.0]))]
xticklabels = [round(xmin, 1), round((xmin + xmax)/2, 1), round(xmax, 1)]
yticklabels = [round(ymin, 1), round((ymin + ymax)/2, 1), round(ymax, 1)]
extent = 0, nc, 0, nr

zeros = np.zeros([nr, nc])
a20 = np.full([nr, nc], a2low)
a21 = np.full([nr, nc], a2high)

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 18))
fig.delaxes(axs[0, 1])
fig.supxlabel(xlabel, x=0.515, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.04, y=0.510, fontsize=fslabel)

letter = ord('b')
for axrow in axs:
    for ax in axrow:
        ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
        if ax.get_subplotspec().is_first_row():
            ax.set_title('Game types', pad=50.0, fontsize=fslabel)
            ax.text(0, nr*1.035, 'a', fontsize=fslabel, weight='bold')
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        else:
            ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
            letter += 1
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

for rho in rhos:

    RR = np.full([nr, nc], rho)
    w00 = mymodule.fitness(a20, a20, zeros, AA, RR)
    w11 = mymodule.fitness(a21, a21, zeros, AA, RR)
    a2social = np.copy(a20)
    mask = (w11 > w00)
    a2social[mask] = a21[mask]
    wsocial = mymodule.fitness(a2social, a2social, zeros, AA, RR)
    Z = np.full([nr, nc, 4], mymodule.colormap['default'])
    a2eq = np.copy(zeros)
    weq = np.copy(zeros)
    xeq = np.copy(zeros)
    w01 = mymodule.fitness(a20, a21, GG, AA, RR)
    w10 = mymodule.fitness(a21, a20, GG, AA, RR)

    mask0 = (w00 > w11)
    T = np.copy(w01)
    R = np.copy(w00)
    P = np.copy(w11)
    S = np.copy(w10)
    mymodule.gametypes(mask0, T, R, P, S, a20, a21, Z, a2eq, xeq, weq)

    mask0 = (w00 < w11)
    T = np.copy(w10)
    R = np.copy(w11)
    P = np.copy(w00)
    S = np.copy(w01)
    mymodule.gametypes(mask0, T, R, P, S, a21, a20, Z, a2eq, xeq, weq)

    Mss = [[a2social, wsocial], [a2eq, weq]]

    axs[0, 0].imshow(Z, extent=extent)

    for row, Ms in zip(axs[1:], Mss):
        for ax, M, traitvmax in zip(row, Ms, traitvmaxs):
            ax.imshow(M, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)

    if movie:
        loges = -log(1.0 - rho, 2)
        text = fig.text(0.90, 0.90, f'loges\n{loges:4.2f}', fontsize=fstick, color='grey', ha='right')
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
