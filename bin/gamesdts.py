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

traitvmaxs = [mymodule.a2max, mymodule.fitness(np.array([mymodule.a2max]), np.array([mymodule.a2max]), np.array([0.0]), np.array([0.9]), np.array([5.0]))]
xticklabels = [round(xmin), round((xmin + xmax)/2), round(xmax)]
yticklabels = [round(ymin, 1), round((ymin + ymax)/2, 1), round(ymax, 1)]
extent = 0, nc, 0, nr

zeros = np.zeros([nr, nc])
a20 = np.copy(zeros)
a21 = a20 + mymodule.a2max/2.0
a22 = a20 + mymodule.a2max
w00 = mymodule.fitness(a20, a20, zeros, AA, RR)
w11 = mymodule.fitness(a21, a21, zeros, AA, RR)
w22 = mymodule.fitness(a22, a22, zeros, AA, RR)
a2social = np.copy(zeros)
mask = (w11 > w00)
a2social[mask] = a21[mask]
mask = (w22 > w11)
a2social[mask] = a22[mask]
wsocial = mymodule.fitness(a2social, a2social, zeros, AA, RR)

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 18))
fig.supxlabel(xlabel, x=0.515, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.04, y=0.510, fontsize=fslabel)

letter = ord('a')
for axrow in axs:
    for ax in axrow:
        ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
        letter += 1
        ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
        if ax.get_subplotspec().is_first_row():
            ax.set_title('Game types', pad=50.0, fontsize=fslabel)
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

    a2eq0 = np.copy(zeros)
    a2eq1 = np.copy(zeros)
    weq0 = np.copy(zeros)
    weq1 = np.copy(zeros)
    xeq = np.copy(zeros)
    w01 = mymodule.fitness(a20, a21, given, AA, RR)
    w10 = mymodule.fitness(a21, a20, given, AA, RR)
    w12 = mymodule.fitness(a21, a22, given, AA, RR)
    w21 = mymodule.fitness(a22, a21, given, AA, RR)

    Z = np.full([nr, nc, 4], mymodule.colormap['default'])
    TS = np.zeros([nr, nc])

    mask0 = (w00 > w11)
    T = np.copy(w01)
    R = np.copy(w00)
    P = np.copy(w11)
    S = np.copy(w10)
    mymodule.gametypes(mask0, T, R, P, S, a20, a21, Z, TS, a2eq0, xeq, weq0)

    mask0 = (w00 < w11)
    T = np.copy(w10)
    R = np.copy(w11)
    P = np.copy(w00)
    S = np.copy(w01)
    mymodule.gametypes(mask0, T, R, P, S, a21, a20, Z, TS, a2eq0, xeq, weq0)

    axs[0, 0].imshow(TS, extent=extent, cmap='magma', vmin=0, vmax=0.7)

    Z = np.full([nr, nc, 4], mymodule.colormap['default'])
    TS = np.zeros([nr, nc])

    mask0 = (w11 > w22)
    T = np.copy(w12)
    R = np.copy(w11)
    P = np.copy(w22)
    S = np.copy(w21)
    mymodule.gametypes(mask0, T, R, P, S, a21, a22, Z, TS, a2eq1, xeq, weq1)

    mask0 = (w11 < w22)
    T = np.copy(w21)
    R = np.copy(w22)
    P = np.copy(w11)
    S = np.copy(w12)
    mymodule.gametypes(mask0, T, R, P, S, a22, a21, Z, TS, a2eq1, xeq, weq1)

    axs[0, 1].imshow(TS, extent=extent, cmap='magma', vmin=0, vmax=0.7)

    a2eq = np.copy(a2eq0)
    weq = np.copy(weq0)
    mask = (a2eq0 == a21)
    a2eq[mask] = a2eq1[mask]
    weq[mask] = weq1[mask]

    Mss = [[a2social, wsocial], [a2eq, weq]]

    for row, Ms in zip(axs[1:], Mss):
        for ax, M, traitvmax in zip(row, Ms, traitvmaxs):
            ax.imshow(M, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)

    if movie:
        text = fig.text(0.90, 0.90, f'given: {given:4.2f}', fontsize=fstick, color='grey', ha='right')
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