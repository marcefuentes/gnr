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
givenmin = 0.0
givenmax = 1.0

num = 21    # Number of subplot rows and columns
ngiven = 21
filename = 'gamesc'

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def gametype (x, y, Z, given, alpha, rho):
    T = mymodule.fitness(x, y, given, AA, RR)
    R = mymodule.fitness(x, x, given, AA, RR)
    P = mymodule.fitness(y, y, given, AA, RR)
    S = mymodule.fitness(y, x, given, AA, RR)
    mask = (R < P)
    H = T[mask]
    T[mask] = S[mask]
    S[mask] = H
    H = R[mask]
    R[mask] = P[mask]
    P[mask] = H
    Z[(T < R) & (P < S)] = mymodule.nodilemma
    Z[(T < R) & (P < S) & (2.0*R <= T + S)] = mymodule.RTSnd
    Z[(T > R) & (P > S)] = mymodule.prisoner
    Z[(T > R) & (P > S) & (2.0*R <= T + S)] = mymodule.RTSpd
    Z[(T >= R) & (P <= S)] = mymodule.snowdrift
    Z[(T >= R) & (P <= S) & (2.0*R <= T + S)] = mymodule.RTSsd
    Z[R == P] = mymodule.nodilemma
    return Z

if givenmin != givenmax:
    movie = True
    givens = np.linspace(givenmin, givenmax, num=ngiven)
    frames = []
else:
    movie = False 
    givens = np.array([givenmin])

nc = num
nr = num
MRT0 = mymodule.b*mymodule.Rq
if givens[-1] > 0.9999999:
    givens[-1] = 0.9999999
alphas = np.linspace(alphamax, alphamin, num=nr)
logess = np.linspace(logesmin, logesmax, num=nc)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
zeros = np.full([nr, nc], 0.0)
ones = np.full([nr, nc], 1.0)
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
extenta2 = 0, nc, 0, nr*2

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 18))
fig.delaxes(axs[0, 1])
fig.supxlabel(xlabel, x=0.515, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.04, y=0.510, fontsize=fslabel)

letter = ord('b')
for axrow in axs:
    for ax in axrow:
        if ax.get_subplotspec().is_first_row():
            ax.set(xticks=[0, nc/2, nc], yticks=[0, nr, nr*2], xticklabels=[], yticklabels=[])
            ax.set_title('Game types', pad=50.0, fontsize=fslabel)
            ax.text(0, nr*2*1.035, 'a', fontsize=fslabel, weight='bold')
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

for given in givens:

    if movie:
        text = fig.text(0.90, 0.90, f'given: {given:4.2f}', fontsize=fstick, color='grey', ha='right')

    MRT = MRT0*(1.0 - given)
    Q0 = mymodule.Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)
    Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
    wsocial = mymodule.fitness(a2social, a2social, given, AA, RR)
    weq = mymodule.fitness(a2eq, a2eq, given, AA, RR)
    Mss = [[a2social, wsocial], [a2eq, weq]]

    Z0 = np.full([nr, nc, 4], mymodule.colormap['default'])
    a2 = 1.01*a2eq
    mask = (a2 > 1.0)
    a2[mask] = ones[mask]
    Z0 = gametype(a2eq, a2, Z0, given, AA, RR)

    Z1 = np.full([nr, nc, 4], mymodule.colormap['default'])
    a2 = 0.99*a2eq
    mask = (a2 < 0.0)
    a2[mask] = zeros[mask]
    Z1 = gametype(a2eq, a2, Z1, given, AA, RR)

    Z = np.full([nr*2, nc, 4], mymodule.colormap['default'])
    Z[::2,:] = Z1
    Z[1::2,:] = Z0

    axs[0, 0].imshow(Z, extent=extenta2, aspect=0.5)

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
