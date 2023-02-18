#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter ()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

traitlabels = ['Games', 'Effort to get $\it{B}$', 'Fitness']
givens = [0.95, 0.50]

num = 21    # Number of subplot rows and columns

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=num)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
MRT0 = mymodule.b*mymodule.Rq
Q0 = mymodule.Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)

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

fig, axs = plt.subplots(nrows=len(givens),
                        ncols=len(traitlabels),
                        figsize=(6*len(traitlabels), 6*(len(givens))))
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

for axrow, given in zip(axs, givens):

    MRT = MRT0*(1.0 - given)
    Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
    weq = mymodule.fitness(a2eq, a2eq, given, AA, RR)
    low = np.full([num, num], a2eq)
    high = np.full([num, num], a2social)
    mask = (a2eq > a2social)
    low[mask] = a2social[mask]
    high[mask] = a2eq[mask]
    T = mymodule.fitness(high, low, given, AA, RR)
    S = mymodule.fitness(low, high, given, AA, RR)
    R = mymodule.fitness(high, high, given, AA, RR)
    P = mymodule.fitness(low, low, given, AA, RR)

    Z = np.full([num, num, 4], mymodule.colormap['default'])
    mymodule.gametypes(T, R, P, S, Z)
    axrow[0].imshow(Z, extent=extent)
    #axrow[1].imshow(a2eq, extent=extent, cmap='viridis', vmin=0, vmax=traitvmaxs[0])
    #axrow[2].imshow(weq, extent=extent, cmap='viridis', vmin=0, vmax=traitvmaxs[1])
    axrow[1].imshow(low, extent=extent, cmap='viridis', vmin=0, vmax=1)
    axrow[2].imshow(high, extent=extent, cmap='viridis', vmin=0, vmax=1)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
