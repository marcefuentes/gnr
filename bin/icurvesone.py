#! /usr/bin/env python

from matplotlib import cm
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

givens = np.linspace(0.0, 1.0, num=21)

num = 3    # Number of subplot rows & columns
numa2 = 128
n_ic = 5    # Number of indifference curves

plotsize = 6

def indifference(q, w, alpha, rho):
    if rho == 0.0:
        if q == 0.0:
            q2 = 1000.0
        else:
            q2 = pow(w/pow(q, 1.0 - alpha), 1.0/alpha)
    elif rho < 0.0:
        if q == 0.0:
            q2 = 1000.0
        else:
            if pow(w, rho) <= (1.0 - alpha)*pow(q, rho):
                q2 = 1000.0
            else:
                q2 = pow((pow(w, rho) - (1.0 - alpha)*pow(q, rho))/
                        alpha, 1.0/rho)
    else:
        if pow(w, rho) <= (1.0 - alpha)*pow(q, rho):
            q2 = -0.1
        else:
            q2 = pow((pow(w, rho) - (1.0 - alpha)*pow(q, rho))/
                    alpha, 1.0/rho)
    return q2

alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=num)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
a1_budget = np.linspace(0.0, mymodule.a1max, num=3)
q2_budget = (mymodule.a2max - mymodule.b*a1_budget)*mymodule.R2
q1_budget = a1_budget*mymodule.R1
q1_ic = np.linspace(0.0, mymodule.a1max*mymodule.R1, num=numa2)
RR, AA = np.meshgrid(rhos, alphas)
ws = np.linspace(2.0/(n_ic + 1), 2.0*n_ic/(n_ic + 1), num=n_ic)
icsss = []
for alpha in alphas:
    icss = []
    for rho in rhos:
        ics = []
        for w in ws:
            ic = np.zeros(numa2)
            for i, q1 in enumerate(q1_ic):
                ic[i] = indifference(q1, w, alpha, rho)
            ics.append(ic)
        icss.append(ics)
    icsss.append(icss)

xlim=[0.0, mymodule.a1max*mymodule.R1]
ylim=[0.0, mymodule.a2max*mymodule.R2]
step = int(num/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
traitvmax = mymodule.fitness(np.array([mymodule.a2max]),
                             np.array([mymodule.a2max]),
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

frames = []
for given in givens:

    a2eq = mymodule.a2eq(given, AA, RR)
    w = mymodule.fitness(a2eq, a2eq, given, AA, RR)
    q2 = a2eq*mymodule.R2
    q2b = q2_budget*(1.0 - given)

    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            for line in axs[i, j].get_lines():
                line.remove()
            for n in range(n_ic): 
                axs[i, j].plot(q1_ic, icsss[i][j][n], c='0.850')
            budget = q2b + q2[i, j]*given
            axs[i, j].plot(q1_budget, budget, c='black', alpha=0.8)
            y = []
            for q1 in q1_ic:
                y.append(indifference(q1, w[i, j], alpha, rho))
            axs[i, j].plot(q1_ic,
                           y,
                           linewidth=4,
                           alpha=0.8,
                           c=cm.viridis(w[i, j]/traitvmax))
    text = fig.text(0.90,
                    0.90,
                    'Given: ' + f'{given:4.2f}',
                    fontsize=biglabels,
                    color='grey',
                    ha='right')
    plt.savefig('temp.png', transparent=False)
    text.remove()
    frames.append(iio.imread('temp.png'))
    os.remove('temp.png')

plt.close()

iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
