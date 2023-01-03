#! /usr/bin/env python

from matplotlib import cm
import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import time

start_time = time.perf_counter ()

alphamin = 0.1
alphamax = 0.9
logesmin = -5.0
logesmax = 5.0
givenmin = 0.95
givenmax = 0.95

num = 21    # Number of subplot rows and columns
numa2 = 3
ngiven = 21
filename = 'landscaped'

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

every = int(num/2)

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
a2 = np.linspace(0.0, mymodule.a2max, num=numa2)
RR, AA = np.meshgrid(rhos, alphas)

xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'

traitvmax = mymodule.fitness(np.array([mymodule.a2max]), np.array([mymodule.a2max]), np.array([0.0]), np.array([0.9]), np.array([5.0]))

zeros = np.zeros([nr, nc])
a20 = np.copy(zeros)
a21 = np.full([nr, nc], mymodule.a2max/2.0)
a22 = np.full([nr, nc], mymodule.a2max)
w00 = mymodule.fitness(a20, a20, zeros, AA, RR)
w11 = mymodule.fitness(a21, a21, zeros, AA, RR)
w22 = mymodule.fitness(a22, a22, zeros, AA, RR)
a2social = np.copy(zeros)
mask = (w11 > w00)
a2social[mask] = a21[mask]
mask = (w22 > w11)
a2social[mask] = a22[mask]
wsocial = mymodule.fitness(a2social, a2social, zeros, AA, RR)

fig = plt.figure(figsize=(8, 8))
fig.supxlabel(xlabel, x=0.56, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.05, y=0.52, fontsize=fslabel)
grid = fig.add_gridspec(nrows=nr, ncols=nc, left=0.22, right=0.9, top=0.86, bottom=0.176, wspace=0, hspace=0)
axs = grid.subplots()

for row in axs:
    for ax in row:
        ax.set(xticks=[], yticks=[], xlim=(0.0, mymodule.a2max), ylim=(0.0, traitvmax))
for ax, loges in zip(axs[-1, ::every], logess[::every]):
    ax.set_xlabel(round(loges), fontsize=fstick)
for ax, alpha in zip(axs[::every, 0], alphas[::every]):
    ax.set_ylabel(f'{alpha:1.1f}', rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

for given in givens:

    if movie:
        text = fig.text(0.90, 0.90, f'given: {given:4.2f}', fontsize=fstick, color='grey', ha='right')

    Z0 = np.full([nr, nc, 4], mymodule.default)
    Z1 = np.full([nr, nc, 4], mymodule.default)
    a2eq0 = np.copy(zeros)
    a2eq1 = np.copy(zeros)
    weq0 = np.copy(zeros)
    weq1 = np.copy(zeros)
    xeq = np.copy(zeros)
    w01 = mymodule.fitness(a20, a21, given, AA, RR)
    w10 = mymodule.fitness(a21, a20, given, AA, RR)
    w12 = mymodule.fitness(a21, a22, given, AA, RR)
    w21 = mymodule.fitness(a22, a21, given, AA, RR)
    w02 = mymodule.fitness(a20, a22, given, AA, RR)
    w20 = mymodule.fitness(a22, a20, given, AA, RR)

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

    for rowax, alpha, rowa2eq, rowweq in zip(axs, alphas, a2eq, weq):
        for ax, rho, cola2eq, colweq in zip(rowax, rhos, rowa2eq, rowweq):
            for line in ax.get_lines():
                line.remove()
            w = mymodule.fitness(np.full(numa2, cola2eq), a2, np.full(numa2, given), np.full(numa2, alpha), np.full(numa2, rho))
            ax.plot(a2, w, linewidth=4, c=cm.magma(colweq/traitvmax))

    if movie:
        plt.savefig('temp.png', transparent=False)
        frames.append(iio.imread('temp.png'))
        os.remove('temp.png')
        text.remove()
    else:
        plt.savefig(filename + '.png', transparent=False)

plt.close()

if movie:
    iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')