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
numa2 = 5
ngiven = 21
filename = 'gamesc'

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
MRT0 = mymodule.b*mymodule.Rq
if givens[-1] > 0.9999999:
    givens[-1] = 0.9999999
alphas = np.linspace(alphamax, alphamin, num=nr)
logess = np.linspace(logesmin, logesmax, num=nc)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
a2 = np.linspace(mymodule.a2max, 0.0, num=numa2)
a2high = np.tile(a2[:-1], [nr, nc])
a2low = np.tile(a2[1:], [nr, nc])
a2high = a2high.T
a2low = a2low.T
RRR, AAA = np.meshgrid(rhos, np.repeat(alphas, numa2-1))
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
extenta2 = 0, nc, 0, nr*(numa2-1)

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 18))
fig.delaxes(axs[0, 1])
fig.supxlabel(xlabel, x=0.515, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.04, y=0.510, fontsize=fslabel)

letter = ord('b')
for axrow in axs:
    for ax in axrow:
        if ax.get_subplotspec().is_first_row():
            ax.set(xticks=[0, nc/2, nc], yticks=[0, nr*(numa2-1)/2, nr*(numa2-1)], xticklabels=[], yticklabels=[])
            ax.set_title('Game types', pad=50.0, fontsize=fslabel)
            ax.text(0, nr*(numa2-1)*1.035, 'a', fontsize=fslabel, weight='bold')
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

    Z = np.full([nr*numa2, nc*numa2, 4], mymodule.colormap['default'])
    T = mymodule.fitness(a2high, a2low, given, AAA, RRR)
    R = mymodule.fitness(a2high, a2high, given, AAA, RRR)
    P = mymodule.fitness(a2low, a2low, given, AAA, RRR)
    S = mymodule.fitness(a2low, a2high, given, AAA, RRR)
    mask = (R < P)
    H = T[mask]
    T[mask] = S[mask]
    S[mask] = H
    H = R[mask]
    R[mask] = P[mask]
    P[mask] = H
    Z[(T < R) & (P < S)] = mymodule.colormap['nodilemma']
    Z[(T < R) & (P < S) & (2.0*R <= T + S)] = mymodule.colormap['nodilemmaRS']
    Z[(T > R) & (P > S)] = mymodule.colormap['prisoner']
    Z[(T > R) & (P > S) & (2.0*R <= T + S)] = mymodule.colormap['prisonerRS']
    Z[(T >= R) & (P <= S)] = mymodule.colormap['snowdrift']
    Z[(T >= R) & (P <= S) & (2.0*R <= T + S)] = mymodule.colormap['snowdriftRS']
    Z[R == P] = mymodule.colormap['nodilemma']

    MRT = MRT0*(1.0 - given)
    Q0 = mymodule.Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)
    Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
    wsocial = mymodule.fitness(a2social, a2social, given, AA, RR)
    weq = mymodule.fitness(a2eq, a2eq, given, AA, RR)
    Mss = [[a2social, wsocial], [a2eq, weq]]

    axs[0, 0].imshow(Z, extent=extenta2, aspect=1.0/(numa2 - 1))

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
