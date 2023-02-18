#! /usr/bin/env python

import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter ()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

traitlabels = ['Games',
                '$\it{T}$ + $\it{S}$ - 2$\it{R}$',
                'Effort to get $\it{B}$',
                'Fitness']
givens = np.linspace(0.0, 1.0, num=20)
a2lows = [0.50, 0.25, 0.00]
rows = a2lows
title = 'Given: '

num = 1001    # Number of subplot rows and columns

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=num)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)

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
extent2 = 0, num, 0, num

fig, axs = plt.subplots(nrows=len(rows),
                        ncols=len(traitlabels),
                        figsize=(6*len(traitlabels), 6*(len(rows))))
fig.supxlabel(xlabel, x=0.513, y=0.01, fontsize=fslabel*1.2)
fig.supylabel(ylabel, x=0.03, y=0.493, fontsize=fslabel*1.2)

letter = ord('a')
for axrow in axs:
    for ax, traitlabel in zip(axrow, traitlabels):
        ax.text(0, 
                num*1.035,
                chr(letter),
                fontsize=fslabel*0.8,
                weight='bold')
        letter += 1
        ax.set(xticks=[0, num/2, num],
                yticks=[0, num/2, num],
                xticklabels=[],
                yticklabels=[])
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
        if ax.get_subplotspec().is_first_row():
            ax.set_title(traitlabel, pad=40.0, fontsize=fslabel*0.9)
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)

frames = []
for given in givens:

    for axrow, a2low in zip(axs, a2lows):

        low = np.full([num, num], a2low)
        high = low + 0.5
        T = mymodule.fitness(high, low, given, AA, RR)
        R = mymodule.fitness(high, high, given, AA, RR)
        P = mymodule.fitness(low, low, given, AA, RR)
        S = mymodule.fitness(low, high, given, AA, RR)
        Z = np.full([num, num, 4], mymodule.colormap['red'])
        mymodule.gamecolors(T, R, P, S, Z)
        axrow[0].imshow(Z, extent=extent2)

        Z = np.zeros([num, num])
        mask = mymodule.prisoner(T, R, P, S) & (2.0*R > T + S)
        Z[mask] = 1.0 - (2.0*R[mask] - T[mask] - S[mask])
        Z = np.ma.masked_where(Z == 0.0, Z)
        cmap = plt.cm.viridis
        cmap.set_bad(color='white')
        axrow[1].imshow(Z, extent=extent2, cmap=cmap, vmin=0, vmax=traitvmaxs[0])

        a2eq = np.zeros([num, num])
        weq = np.zeros([num, num])
        mymodule.equilibrium(T, R, P, S, low, high, a2eq, weq)
        if given == 0.0:
            mask = (R == P)
            a2eq[mask] = (low[mask] + high[mask])/2.0
            weq[mask] = (mymodule.fitness(low[mask], low[mask], given, AA[mask], RR[mask]) +
                            mymodule.fitness(high[mask], high[mask], given, AA[mask], RR[mask]))/2.0
        axrow[2].imshow(a2eq, extent=extent2, cmap='viridis', vmin=0, vmax=traitvmaxs[0])
        axrow[3].imshow(weq, extent=extent2, cmap='viridis', vmin=0, vmax=traitvmaxs[1])

    movieframe = given
    text = fig.text(0.90,
                    0.02,
                    title + f'{movieframe:4.2f}',
                    fontsize=fslabel,
                    color='grey',
                    ha='right')
    plt.savefig('temp.png', transparent=False)
    text.remove()
    frames.append(iio.imread('temp.png'))
    os.remove('temp.png')

plt.close()

iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
