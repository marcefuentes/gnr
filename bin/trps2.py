#! /usr/bin/env python

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
a2low = 0.0
a2high = 0.49

num = 21    # Number of subplot rows and columns
ngiven = 21
filename = 'trps'

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
RR, AA = np.meshgrid(rhos, alphas)

xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'

zeros = np.zeros([nr, nc])
a20 = np.full([nr, nc], a2low)
a21 = np.full([nr, nc], a2high)
w00 = fitness(a20, a20, zeros, AA, RR)
w11 = fitness(a21, a21, zeros, AA, RR)
a2social = np.copy(a20)
mask = (w11 > w00)
a2social[mask] = a21[mask]
wsocial = fitness(a2social, a2social, zeros, AA, RR)
xaxis = [1, 2, 3, 4]

fig = plt.figure(figsize=(8, 8))
fig.supxlabel(xlabel, x=0.56, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.05, y=0.52, fontsize=fslabel)
grid = fig.add_gridspec(nrows=nr, ncols=nc, left=0.22, right=0.9, top=0.86, bottom=0.176, wspace=0, hspace=0)
axs = grid.subplots()

for row in axs:
    for ax in row:
        ax.set(xticks=[], yticks=[], xlim=(0, 5), ylim=(0.0, 2.0))
for ax, loges in zip(axs[-1, ::every], logess[::every]):
    ax.set_xlabel(round(loges), fontsize=fstick)
for ax, alpha in zip(axs[::every, 0], alphas[::every]):
    ax.set_ylabel(f'{alpha:1.1f}', rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

for given in givens:

    if movie:
        text = fig.text(0.90, 0.90, f'given: {given:4.2f}', fontsize=fstick, color='grey', ha='right')

    Z = np.full([nr, nc, 4], mymodule.default)
    Tf = np.copy(zeros)
    Rf = np.copy(zeros)
    Pf = np.copy(zeros)
    Sf = np.copy(zeros)
    a2eq = np.copy(zeros)
    weq = np.copy(zeros)
    xeq = np.copy(zeros)
    w01 = fitness(a20, a21, given, AA, RR)
    w10 = fitness(a21, a20, given, AA, RR)

    mask0 = (w00 > w11)
    T = np.copy(w01)
    R = np.copy(w00)
    P = np.copy(w11)
    S = np.copy(w10)
    Tf[mask0] = T[mask0]
    Rf[mask0] = R[mask0]
    Pf[mask0] = P[mask0]
    Sf[mask0] = S[mask0]
    mymodule.gametypes(mask0, T, R, P, S, a20, a21, Z, a2eq, xeq, weq)

    mask0 = (w00 < w11)
    T = np.copy(w10)
    R = np.copy(w11)
    P = np.copy(w00)
    S = np.copy(w01)
    Tf[mask0] = T[mask0]
    Rf[mask0] = R[mask0]
    Pf[mask0] = P[mask0]
    Sf[mask0] = S[mask0]
    mymodule.gametypes(mask0, T, R, P, S, a21, a20, Z, a2eq, xeq, weq)

    for row, rowT, rowR, rowP, rowS, rowZ in zip(axs, Tf, Rf, Pf, Sf, Z):
        for ax, tt, rr, pp, ss, zz in zip(row, rowT, rowR, rowP, rowS, rowZ):
            y = [tt, rr, pp, ss]
            for line in ax.get_lines():
                line.remove()
            ax.plot(xaxis, y, c=zz, linewidth=3, marker='o', markerfacecolor='white', markersize=3)

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
