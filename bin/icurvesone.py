#! /usr/bin/env python

import os
import time

from matplotlib import cm
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

givens = np.linspace(0.0, 1.0, num=41)

num = 3     # Number of subplot rows & columns
numa2 = 256 # Number of points along each curve
n_ic = 5    # Number of indifference curves

plotsize = 6

def update(given):
    a2private = my.a2eq(given, AA, RR)
    w = my.fitness(a2private, a2private, given, AA, RR)
    q2 = a2private*my.R2
    q2b = q2_budget*(1.0 - given)

    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            new_q2_budget = q2b + q2[i, j]*given
            budget[i, j].set_ydata(new_q2_budget)
            for k, q1 in enumerate(q1_ic):
                ic[k] = my.indifference(q1, w[i, j], alpha, rho)
            icurve[i, j].set_ydata(ic)
            icurve[i, j].set_color(cm.viridis(w[i, j]/traitvmax))
    if len(givens) > 1:
        axs[0, 2].set_title('Given: ' + f'{given:4.2f}',
                            fontsize=ticklabels,
                            color='grey',
                            ha='right',
                            pad=10)
    return np.concatenate([budget.flatten(), icurve.flatten()])
    
alphas = np.linspace(my.alphamax, my.alphamin, num=num)
logess = np.linspace(my.logesmin, my.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
a1_budget = np.array([0.0, my.a1max])
q1_budget = a1_budget*my.R1
q2_budget = (my.a2max - my.b*a1_budget)*my.R2
q1_ic = np.linspace(0.001*my.R1,
                    (my.a1max - 0.001)*my.R1,
                    num=numa2)

RR, AA = np.meshgrid(rhos, alphas)
ws = np.linspace(2.0/(n_ic + 1), 2.0*n_ic/(n_ic + 1), num=n_ic)
ics = np.empty((num, num, n_ic, numa2), dtype=np.float64)
for i, alpha in enumerate(alphas):
    for j, rho in enumerate(rhos):
        for k, w in enumerate(ws):
            for l, q1 in enumerate(q1_ic):
                ics[i, j, k, l] = my.indifference(q1, w, alpha, rho)

xlim=[0.0, my.a1max*my.R1]
ylim=[0.0, my.a2max*my.R2]
step = int(num/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
traitvmax = my.fitness(np.array([my.a2max]),
                       np.array([my.a2max]),
                       np.array([0.0]),
                       np.array([0.9]),
                       np.array([5.0]))

width = plotsize
height = plotsize
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.56,
              y=0.03,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.05,
              y=0.52,
              fontsize=biglabels)

grid = fig.add_gridspec(nrows=num,
                        ncols=num,
                        left=0.22,
                        right=0.9,
                        top=0.86,
                        bottom=0.176,
                        wspace=0,
                        hspace=0)
axs = grid.subplots()

for i, alpha in enumerate(alphas):
    for j, rho in enumerate(rhos):
        axs[i, j].set(xticks=[],
                      yticks=[],
                      xlim=xlim,
                      ylim=ylim)
for i in range(0, num, step):
    axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                         rotation='horizontal',
                         horizontalalignment='right',
                         verticalalignment='center',
                         fontsize=ticklabels)
for j in range(0, num, step):
    axs[-1, j].set_xlabel(f'{logess[j]:2.0f}', fontsize=ticklabels)

ic = np.empty(q1_ic.shape, dtype=np.float64)
budget = np.empty(axs.shape, dtype=object)
icurve = np.empty(axs.shape, dtype=object)

for i, alpha in enumerate(alphas):
    for j, rho in enumerate(rhos):
        for k in range(n_ic): 
            axs[i, j].plot(q1_ic, ics[i, j, k], c='0.850')
        budget[i, j], = axs[i, j].plot(q1_budget, q2_budget, c='black', alpha=0.8)
        icurve[i, j], = axs[i, j].plot(q1_ic, ic, linewidth=4, alpha=0.8)

ani = FuncAnimation(fig, update, givens, blit=True)

if len(givens) > 1:
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
