#! /usr/bin/env python

import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import numpy as np
import time

start_time = time.perf_counter ()

minalpha = 0.1
maxalpha = 0.9
minloges = -5.0
maxloges = 5.0
mingiven = 0.95
maxgiven = 0.95

num = 21    # Number of subplot rows and columns
filename = 'games_d_es_alpha'
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

traitlabels = ['Effort to get $\it{B}$', 'Fitness']
traitvmaxs = [1.0, 1.8]
fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def fitness(x, y, given, alpha, rho):
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - given) + x*R2*given
    w = q1*q2
    mask = (w > 0.0) & (rho == 0.0)
    w[mask] = pow(q1[mask], 1.0 - alpha[mask])*pow(q2[mask], alpha[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = (1.0 - alpha[mask])*pow(q1[mask], rho[mask]) + alpha[mask]*pow(q2[mask], rho[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = pow(w[mask], 1.0/rho[mask])
    mask = (rho > 0.0)
    w[mask] = pow((1.0 - alpha[mask])*pow(q1[mask], rho[mask]) + alpha[mask]*pow(q2[mask], rho[mask]), 1.0/rho[mask])
    return w

def gametypes(a2c, a2d):
    mask = (mask0 & (T < R) & (P < S))
    Z[mask] = nodilemma
    a2eq[mask] = a2c[mask]
    weq[mask] = R[mask]
    mask = (mask0 & (T >= R) & (P <= S) & (R != P))
    Z[mask] = snowdrift
    xeq[mask] = (P[mask] - S[mask])/(R[mask] - S[mask] - T[mask] + P[mask])
    a2eq[mask] = a2c[mask]*xeq[mask] + a2d[mask]*(1.0 - xeq[mask])
    weq[mask] = (T[mask] + S[mask])*xeq[mask]*(1.0 - xeq[mask]) + R[mask]*xeq[mask]*xeq[mask] + P[mask]*(1.0 - xeq[mask])*(1.0 - xeq[mask])
    mask = (mask0 & (T > R) & (P > S))
    Z[mask] = prisoner
    a2eq[mask] = a2d[mask]
    weq[mask] = P[mask]
    pass

if mingiven != maxgiven:
    movie = True
    givens = np.linspace(mingiven, maxgiven, num=num)
    frames = []
else:
    movie = False 
    givens = np.array([mingiven])

b = a2max/a1max
nc = num
nr = num
alphas = np.linspace(maxalpha, minalpha, num=nr)
logess = np.linspace(minloges, maxloges, num=nc)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)

minx = minloges
maxx = maxloges
xlabel = 'Substitutability of $\it{B}$'
miny = minalpha
maxy = maxalpha
ylabel = 'Value of $\it{B}$'

xticklabels = [round(minx), round((minx + maxx)/2), round(maxx)]
yticklabels = [miny, (miny + maxy)/2, maxy]
extent = 0, nr, 0, nc
prisoner = [0.5, 0.0, 0.0, 1.0]
snowdrift = [0.0, 1.0, 1.0, 1.0]
nodilemma = [1.0, 1.0, 1.0, 1.0]
green = [0.0, 1.0, 0.0, 1.0]

zeros = np.zeros([nr, nc])
a20 = np.copy(zeros)
a21 = np.full([nr, nc], a2max/2.0)
a22 = np.full([nr, nc], a2max)
w00 = fitness(a20, a20, zeros, AA, RR)
w11 = fitness(a21, a21, zeros, AA, RR)
w22 = fitness(a22, a22, zeros, AA, RR)
a2social = np.copy(zeros)
mask = (w11 > w00)
a2social[mask] = a21[mask]
mask = (w22 > w11)
a2social[mask] = a22[mask]
wsocial = fitness(a2social, a2social, zeros, AA, RR)

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

for given in givens:

    if movie:
        text = fig.text(0.80, 0.80, f'given\n{given:4.2f}', fontsize=fstick+4, color='grey', ha='right')

    Z = np.full([nr, nc, 4], green)
    a2eq = np.copy(zeros)
    weq = np.copy(zeros)
    xeq = np.copy(zeros)
    w01 = fitness(a20, a21, given, AA, RR)
    w10 = fitness(a21, a20, given, AA, RR)
    w12 = fitness(a21, a22, given, AA, RR)
    w21 = fitness(a22, a21, given, AA, RR)
    w02 = fitness(a20, a22, given, AA, RR)
    w20 = fitness(a22, a20, given, AA, RR)

    mask0 = (wsocial == w00)
    T = np.copy(w01)
    R = np.copy(w00)
    P = np.copy(w11)
    S = np.copy(w10)
    gametypes(a20, a21)

    mask0 = (wsocial == w11) & (w20 > w22)
    T = np.copy(w10)
    R = np.copy(w11)
    P = np.copy(w00)
    S = np.copy(w01)
    gametypes(a21, a20)

    mask0 = (wsocial == w11) & (w20 <= w22)
    T = np.copy(w12)
    R = np.copy(w11)
    P = np.copy(w22)
    S = np.copy(w21)
    gametypes(a21, a22)

    mask0 = (wsocial == w22) & (w10 < w11)
    T = np.copy(w21)
    R = np.copy(w22)
    P = np.copy(w11)
    S = np.copy(w12)
    gametypes(a22, a21)

    mask0 = (wsocial == w22) & (w10 >= w11)
    T = np.copy(w10)
    R = np.copy(w11)
    P = np.copy(w00)
    S = np.copy(w01)
    gametypes(a21, a20)

    Mss = [[a2social, wsocial], [a2eq, weq]]

    axs[0, 0].imshow(Z, extent=extent)

    for row, Ms in zip(axs[1:], Mss):
        for ax, M, traitvmax in zip(row, Ms, traitvmaxs):
            ax.imshow(M, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)

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
