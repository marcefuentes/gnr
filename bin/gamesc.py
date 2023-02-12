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
numa2 = 64
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

alphas = np.linspace(alphamax, alphamin, num=num)
logess = np.linspace(logesmin, logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
nr = len(alphas)
nc = len(logess)
RR, AA = np.meshgrid(rhos, alphas)
MRT0 = mymodule.b*mymodule.Rq
Q0 = mymodule.Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)
wsocial = mymodule.fitness(a2social, a2social, 0.0, AA, RR)
X, Y = np.meshgrid(np.linspace(0.0, mymodule.a2max, num=numa2),
                    np.linspace(mymodule.a2max, 0.0, num=numa2))
X = np.tile(A=X, reps=[nr, nc])
Y = np.tile(A=Y, reps=[nr, nc])
H = np.copy(Y)
Y[(X > Y)] = X[(X > Y)] 
X[(X > H)] = H[(X > H)] 
RRR, AAA = np.meshgrid(np.repeat(rhos, numa2),
                        np.repeat(alphas, numa2))
R = mymodule.fitness(Y, Y, 0.0, AAA, RRR)
P = mymodule.fitness(X, X, 0.0, AAA, RRR)

xmin = logess[0]
xmax = logess[-1]
xlabel = 'Substitutability of $\it{B}$'
ymin = alphas[-1]
ymax = alphas[0]
ylabel = 'Value of $\it{B}$'

traitvmaxs = [mymodule.a2max,
                mymodule.fitness(np.array([mymodule.a2max]),
                                    np.array([mymodule.a2max]),
                                    np.array([0.0]),
                                    np.array([0.9]),
                                    np.array([5.0]))]
xticklabels = [round(xmin),
                round((xmin + xmax)/2),
                round(xmax)]
yticklabels = [round(ymin, 1),
                round((ymin + ymax)/2, 1),
                round(ymax, 1)]
extent = 0, nc, 0, nr
extenta2 = 0, nc*numa2, 0, nr*numa2

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 18))
fig.supxlabel(xlabel, x=0.515, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.04, y=0.510, fontsize=fslabel)

letter = ord('a')
for axrow in axs:
    for ax in axrow:
        if ax.get_subplotspec().is_first_row():
            ax.text(0,
                    nr*numa2*1.035,
                    chr(letter),
                    fontsize=fslabel,
                    weight='bold')
            ax.set(xticks=[0, nc*numa2/2, nc*numa2],
                    yticks=[0, nr*numa2/2, nr*numa2],
                    xticklabels=[],
                    yticklabels=[])
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        else:
            ax.set(xticks=[0, nc/2, nc],
                            yticks=[0, nr/2, nr],
                            xticklabels=[],
                            yticklabels=[])
            if letter<= ord('z'):
                ax.text(0, 
                        nr*1.035,
                        chr(letter),
                        fontsize=fslabel,
                        weight='bold')
            else:
                ax.text(0,
                        nr*1.035,
                        'a' + chr(letter - 26),
                        fontsize=fslabel,
                        weight='bold')
            if ax.get_subplotspec().is_first_col():
                ax.set_yticklabels(yticklabels, fontsize=fstick) 
        letter += 1
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

if givens[-1] > 0.9999999:
    givens[-1] = 0.9999999

for given in givens:

    T = mymodule.fitness(Y, X, given, AAA, RRR)
    S = mymodule.fitness(X, Y, given, AAA, RRR)
    Z = np.full([nr*numa2, nc*numa2, 4], mymodule.colormap['default'])
    mymodule.gametypes(T, R, P, S, Z)
    ax = axs[0, 0]
    ax.imshow(Z, extent=extenta2)
    ax.set_title('Game types', pad=50.0, fontsize=fslabel)

    Z = np.zeros([nr*numa2, nc*numa2]) + 1.0
    mask = ((T > R) & (R >= P) & (P > S)) | ((T >= P) & (P > R) & (R >= S) & (2.0*P < T + S))
    Z[mask] = - T[mask] - S[mask] + 2.0*R[mask]
    mask = (((T >= P) & (P > R) & (R >= S) & (2.0*P > T + S))) | ((T < R) & (R > P) & (P < S)) 
    Z = np.ma.masked_array(Z, mask)
    cmap = plt.cm.cool
    cmap.set_bad('white')
    ax = axs[0, 1]
    ax.imshow(Z, extent=extenta2, cmap=cmap, vmin=0.0, vmax=1.0)
    ax.set_title('Value of taking turns',
                    pad=50.0,
                    fontsize=fslabel)

    MRT = MRT0*(1.0 - given)
    Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
    weq = mymodule.fitness(a2eq, a2eq, given, AA, RR)
    Mss = [[a2social, wsocial], [a2eq, weq]]

    for row, Ms in zip(axs[1:], Mss):
        for ax, M, traitvmax in zip(row, Ms, traitvmaxs):
            ax.imshow(M, extent=extent, cmap='viridis', vmin=0, vmax=traitvmax)

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
