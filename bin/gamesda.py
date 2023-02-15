#! /usr/bin/env python

import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import time

start_time = time.perf_counter ()

traitlabels = ['Games', 'Effort to get $\it{B}$', 'Fitness']
given = 0.95

num = 1001    # Number of subplot rows and columns
filename = 'gamesda'

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=num)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
a2lows = np.linspace(0.0, 0.5, num=20)
a2highs = a2lows + 0.5

xmin = mymodule.logesmin
xmax = mymodule.logesmax
xlabel = 'Substitutability of $\it{B}$'
ymin = mymodule.alphamin
ymax = mymodule.alphamax
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
extent = 0, num, 0, num

fig, axs = plt.subplots(nrows=1,
                        ncols=len(traitlabels),
                        figsize=(6*len(traitlabels), 6))
fig.supxlabel(xlabel, x=0.513, y=0.01, fontsize=fslabel*1.2)
fig.supylabel(ylabel, x=0.03, y=0.493, fontsize=fslabel*1.2)

letter = ord('a')
for ax, traitlabel in zip(axs, traitlabels):
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
for a2low, a2high in zip(a2lows, a2highs):

    low = np.full([num, num], a2low)
    high = np.full([num, num], a2high)
    T = mymodule.fitness(high, low, given, AA, RR)
    S = mymodule.fitness(low, high, given, AA, RR)
    R = mymodule.fitness(high, high, given, AA, RR)
    P = mymodule.fitness(low, low, given, AA, RR)

    Z = np.full([num, num, 4], mymodule.colormap['white'])
    mymodule.gamecolors(T, R, P, S, Z)
    axs[0].imshow(Z, extent=extent)

    a2eq = np.zeros([num, num])
    weq = np.zeros([num, num])
    mymodule.equilibrium(T, R, P, S, low, high, a2eq, weq)

    axs[1].imshow(a2eq, extent=extent, cmap='viridis', vmin=0, vmax=traitvmaxs[0])
    axs[2].imshow(weq, extent=extent, cmap='viridis', vmin=0, vmax=traitvmaxs[1])

    text = fig.text(0.90, 0.90, f'{a2high:4.2f}\n{a2low:4.2f}', fontsize=fstick, color='grey', ha='right')
    plt.savefig('temp.png', transparent=False)
    text.remove()
    frames.append(iio.imread('temp.png'))
    os.remove('temp.png')

plt.close()

iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
