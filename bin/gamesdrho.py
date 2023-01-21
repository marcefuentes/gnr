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
logesmin = 0.0
logesmax = 20.0
givenmin = 0.0
givenmax = 1.0

num = 21    # Number of subplot rows and columns
numa2 = 2
nloges = 21
filename = 'gamesdrho'

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
extenta2 = 0, nc, 0, nr*numa2

zeros = np.zeros([nr, nc])
a20 = np.copy(zeros)
a21 = np.full([nr, nc], mymodule.a2max/2.0)
a22 = np.full([nr, nc], mymodule.a2max)

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 18))
fig.delaxes(axs[0, 1])
fig.supxlabel(xlabel, x=0.515, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.04, y=0.510, fontsize=fslabel)

letter = ord('b')
for axrow in axs:
    for ax in axrow:
        if ax.get_subplotspec().is_first_row():
            ax.set(xticks=[0, nc/2, nc], yticks=[0, nr*numa2/2, nr*numa2], xticklabels=[], yticklabels=[])
            ax.set_title('Game types', pad=50.0, fontsize=fslabel)
            ax.text(0, nr*numa2*1.035, 'a', fontsize=fslabel, weight='bold')
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        else:
            ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
            ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
            letter += 1
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

for rho in rhos:

    RR = np.full([nr, nc], rho)
    w00 = mymodule.fitness(a20, a20, zeros, AA, RR)
    w11 = mymodule.fitness(a21, a21, zeros, AA, RR)
    w22 = mymodule.fitness(a22, a22, zeros, AA, RR)
    a2social = np.copy(zeros)
    mask = (w11 > w00)
    a2social[mask] = a21[mask]
    mask = (w22 > w11)
    a2social[mask] = a22[mask]
    wsocial = mymodule.fitness(a2social, a2social, zeros, AA, RR)
    Z0 = np.full([nr, nc, 4], mymodule.colormap['default'])
    Z1 = np.full([nr, nc, 4], mymodule.colormap['default'])
    a2eq0 = np.copy(zeros)
    a2eq1 = np.copy(zeros)
    weq0 = np.copy(zeros)
    weq1 = np.copy(zeros)
    xeq = np.copy(zeros)
    w01 = mymodule.fitness(a20, a21, GG, AA, RR)
    w10 = mymodule.fitness(a21, a20, GG, AA, RR)
    w12 = mymodule.fitness(a21, a22, GG, AA, RR)
    w21 = mymodule.fitness(a22, a21, GG, AA, RR)
    w02 = mymodule.fitness(a20, a22, GG, AA, RR)
    w20 = mymodule.fitness(a22, a20, GG, AA, RR)

    mask0 = (w00 > w11)
    T = np.copy(w01)
    R = np.copy(w00)
    P = np.copy(w11)
    S = np.copy(w10)
    mymodule.gametypes(mask0, T, R, P, S, a20, a21, Z0, a2eq0, xeq, weq0)

    mask0 = (w00 < w11)
    T = np.copy(w10)
    R = np.copy(w11)
    P = np.copy(w00)
    S = np.copy(w01)
    mymodule.gametypes(mask0, T, R, P, S, a21, a20, Z0, a2eq0, xeq, weq0)

    mask0 = (w11 > w22)
    T = np.copy(w12)
    R = np.copy(w11)
    P = np.copy(w22)
    S = np.copy(w21)
    mymodule.gametypes(mask0, T, R, P, S, a21, a22, Z1, a2eq1, xeq, weq1)

    mask0 = (w11 < w22)
    T = np.copy(w21)
    R = np.copy(w22)
    P = np.copy(w11)
    S = np.copy(w12)
    mymodule.gametypes(mask0, T, R, P, S, a22, a21, Z1, a2eq1, xeq, weq1)

    a2eq = np.copy(a2eq0)
    weq = np.copy(weq0)
    mask = (a2eq0 == a21)
    a2eq[mask] = a2eq1[mask]
    weq[mask] = weq1[mask]

    Mss = [[a2social, wsocial], [a2eq, weq]]

    Z = np.full([nr*numa2, nc, 4], mymodule.colormap['default'])
    Z[::2,:] = Z1
    Z[1::2,:] = Z0

    axs[0, 0].imshow(Z, extent=extenta2, aspect=1.0/numa2)

    for row, Ms in zip(axs[1:], Mss):
        for ax, M, traitvmax in zip(row, Ms, traitvmaxs):
            ax.imshow(M, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)

    if movie:
        loges = -log(1.0 - rho, 2)
        text = fig.text(0.90, 0.90, f'loges: {loges:4.2f}', fontsize=fstick, color='grey', ha='right')
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
