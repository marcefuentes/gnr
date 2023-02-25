#! /usr/bin/env python

import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter ()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

titles = ['Games',
            '$\it{R}$ - $\it{P}$',
            '$\it{T}$ + $\it{S}$ - 2$\it{R}$']
givens = np.linspace(0.95, 1.0, num=1)

ext = 1024

fslarge = 32 # Label font size
fssmall = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

if givens[-1] > 0.9999999:
    givens[-1] = 0.9999999
alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=ext)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=ext)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
MRT0 = mymodule.b*mymodule.Rq
Q0 = mymodule.Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)
highs = [np.full([ext, ext], mymodule.a2max), 
        np.full([ext, ext], a2social)] 

rows = highs
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
letterposition = ext*1.035
xmin = logess[0]
xmax = logess[-1]
ymin = alphas[-1]
ymax = alphas[0]
xticks = [0, ext/2, ext]
yticks = [0, ext/2, ext]
xticklabels = [f'{xmin:2.0f}',
                f'{(xmin + xmax)/2.0:2.0f}',
                f'{xmax:2.0f}']
yticklabels = [f'{ymin:3.1f}',
                f'{(ymin + ymax)/2.0:3.1f}',
                f'{ymax:3.1f}']
extent = 0, ext, 0, ext

fig, axs = plt.subplots(nrows=len(rows),
                        ncols=len(titles),
                        figsize=(6*len(titles), 6*len(rows)))
fig.supxlabel(xlabel,
                x=0.513,
                y=0.01,
                fontsize=fslarge*1.2)
fig.supylabel(ylabel,
                x=0.03,
                y=0.493,
                fontsize=fslarge*1.2)

for ax in fig.get_axes():
    ax.set(xticks=xticks,
            yticks=yticks,
            xticklabels=[],
            yticklabels=[])
    ax.text(0,
            letterposition,
            chr(letter),
            fontsize=fslarge*0.8,
            weight='bold')
    letter += 1 
for i, row in enumerate(rows):
    axs[i, 0].set_yticklabels(yticklabels, fontsize=fssmall)
for j, title in enumerate(titles):
    axs[0, j].set_title(title, pad=40.0, fontsize=fslarge)
    axs[-1, j].set_xticklabels(xticklabels, fontsize=fssmall)

frames = []
for given in givens:

    MRT = MRT0*(1.0 - given)
    Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
    lows = [np.full([ext, ext], a2eq),
            np.full([ext, ext], 0.0)]

    for i, (low, high) in enumerate(zip(lows, highs)):

        T = mymodule.fitness(high, low, given, AA, RR)
        R = mymodule.fitness(high, high, given, AA, RR)
        P = mymodule.fitness(low, low, given, AA, RR)
        S = mymodule.fitness(low, high, given, AA, RR)
        Z = np.full([ext, ext, 4], mymodule.colormap['red'])

        mymodule.gamecolors(T, R, P, S, Z)
        axs[i, 0].imshow(Z, extent=extent)

        N = np.full([ext, ext, 4], mymodule.colormap['transparent'])
        mymodule.nodilemmacolors(T, R, P, S, N)

        Z = np.zeros([ext, ext])
        Z = R - P
        axs[i, 1].imshow(Z, extent=extent, vmin=-1, vmax=1)
        axs[i, 1].imshow(N, extent=extent)

        Z = np.zeros([ext, ext])
        Z = T + S - 2.0*R
        axs[i, 2].imshow(Z, extent=extent, vmin=-1, vmax=1)
        axs[i, 2].imshow(N, extent=extent)

    text = fig.text(0.90,
                    0.02,
                    'Given: ' + f'{given:4.2f}',
                    fontsize=fslarge,
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
