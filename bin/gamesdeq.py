#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter ()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

traitlabels = ['Games (lower)',
                'Games (upper)',
                '$\it{S}$ - $\it{P}$',
                '$\it{T}$ + $\it{S}$ - 2$\it{R}$']
given = 0.95

num = 1001    # Number of subplot rows and columns

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

if given > 0.9999999:
    given = 0.9999999
alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=num)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
MRT0 = mymodule.b*mymodule.Rq
Q0 = mymodule.Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)
MRT = MRT0*(1.0 - given)
Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
a2lows = [0.0, a2eq]
a2highs = [a2social, mymodule.a2max]

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
    ax.set_title(traitlabel, pad=40.0, fontsize=fslabel*0.9)
    ax.set_xticklabels(xticklabels, fontsize=fstick)

for i, (a2low, a2high) in enumerate(zip(a2lows, a2highs)):

    low = np.full([num, num], a2low)
    high = np.full([num, num], a2high)
    T = mymodule.fitness(high, low, given, AA, RR)
    R = mymodule.fitness(high, high, given, AA, RR)
    P = mymodule.fitness(low, low, given, AA, RR)
    S = mymodule.fitness(low, high, given, AA, RR)
    Z = np.full([num, num, 4], mymodule.colormap['red'])
    mymodule.gamecolors(T, R, P, S, Z)
    axs[2*i].imshow(Z, extent=extent2)

    if i == 0:
        Z = np.zeros([num, num])
        mask = mymodule.prisoner(T, R, P, S)
        Z[mask] = 1.0 - (P[mask] - S[mask])
        Z = np.ma.masked_where(Z == 0.0, Z)
        cmap = plt.cm.viridis
        cmap.set_bad(color='white')
        axs[1].imshow(Z, extent=extent2, cmap=cmap)
    else:
        Z = np.zeros([num, num])
        mask = mymodule.prisoner(T, R, P, S) & (2.0*R > T + S)
        Z[mask] = 1.0 - (2.0*R[mask] - T[mask] - S[mask])
        Z = np.ma.masked_where(Z == 0.0, Z)
        cmap = plt.cm.viridis
        cmap.set_bad(color='white')
        axs[3].imshow(Z, extent=extent2, cmap=cmap)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
