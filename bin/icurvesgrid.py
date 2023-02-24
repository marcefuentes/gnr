#! /usr/bin/env python

from matplotlib import cm
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter ()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

givens = [0.0, 0.5, 0.95]
num = 3    # Number of subplot rows & columns
numa2 = 128
n_ic = 5    # Number of indifference curves

fslarge = 32 # Label font size
fssmall = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

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

if givens[-1] > 0.9999999:
    givens[-1] = 0.9999999
alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=num)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
a1_budget = np.linspace(0.0, mymodule.a1max, num=3)
q2_budget = (mymodule.a2max - mymodule.b*a1_budget)*mymodule.R2
q1_budget = a1_budget*mymodule.R1
q1_ic = np.linspace(0.0, mymodule.a1max*mymodule.R1, num=numa2)
RR, AA = np.meshgrid(rhos, alphas)
MRT0 = mymodule.b*mymodule.Rq
wis = np.linspace(2.0/(n_ic + 1), 2.0*n_ic/(n_ic + 1), num=n_ic)
icsss = []
for alpha in alphas:
    icss = []
    for rho in rhos:
        ics = []
        for w in wis:
            ic = []
            for q1 in q1_ic:
                ic.append(indifference(q1, w, alpha, rho))
            ics.append(ic)
        icss.append(ics)
    icsss.append(icss)

xlim=[0.0, mymodule.a1max*mymodule.R1]
ylim=[0.0, mymodule.a2max*mymodule.R2]
step = int(num/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
traitvmax = mymodule.fitness(np.array([mymodule.a2max]),
                                np.array([mymodule.a2max]),
                                np.array([0.0]),
                                np.array([0.9]),
                                np.array([5.0]))

fig = plt.figure(figsize=(6*len(givens), 6))
fig.supxlabel(xlabel,
                x=0.525,
                y=0.0,
                fontsize=fslarge)
fig.supylabel(ylabel,
                x=0.08,
                y=0.52,
                fontsize=fslarge)

outergrid = fig.add_gridspec(nrows=1,
                                ncols=len(givens),
                                left=0.15,
                                right=0.9,
                                top=0.86,
                                bottom=0.176)

for g, given in enumerate(givens):
    grid = outergrid[g].subgridspec(nrows=num,
                                    ncols=num,
                                    wspace=0,
                                    hspace=0)
    axs = grid.subplots()
    axs[0, int(num/2)].set_title(f'{given*100:2.0f}%',
                                    pad=30.0,
                                    fontsize=fslarge*0.8)
    axs[0, 0].set_title(chr(letter),
                        fontsize=fslarge*0.8,
                        weight='bold',
                        loc='left')
    letter += 1

    for ax in fig.get_axes():
        ax.set(xticks=[],
                yticks=[],
                xlim=xlim,
                ylim=ylim)
    if g == 0:
        for i in range(0, num, step):
            axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                                    rotation='horizontal',
                                    horizontalalignment='right',
                                    verticalalignment='center',
                                    fontsize=fssmall)
    for j in range(0, num, step):
        axs[-1, j].set_xlabel(f'{logess[j]:2.0f}', fontsize=fssmall)

    MRT = MRT0*(1.0 - given)
    Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
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

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
