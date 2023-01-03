#! /usr/bin/env python

from matplotlib import cm
import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import time

start_time = time.perf_counter ()

alphamin = 0.1
alphamax = 0.9
logesmin = -5.0
logesmax = 5.0
givenmin = 0.95
givenmax = 0.95

num = 21    # Number of subplot rows and columns
numa2 = 64
ngiven = 21
filename = 'landscapec'

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

every = int(num/2)

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
a2 = np.linspace(0.0, mymodule.a2max, num=numa2)

xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'

traitvmax = mymodule.fitness(np.array([mymodule.a2max]), np.array([mymodule.a2max]), np.array([0.0]), np.array([0.9]), np.array([5.0]))

fig = plt.figure(figsize=(8, 8))
fig.supxlabel(xlabel, x=0.56, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.05, y=0.52, fontsize=fslabel)
grid = fig.add_gridspec(nrows=nr, ncols=nc, left=0.22, right=0.9, top=0.86, bottom=0.176, wspace=0, hspace=0)
axs = grid.subplots()

for row in axs:
    for ax in row:
        ax.set(xticks=[], yticks=[], xlim=(0.0, mymodule.a2max), ylim=(0.0, traitvmax))
for ax, loges in zip(axs[-1, ::every], logess[::every]):
    ax.set_xlabel(round(loges), fontsize=fstick)
for ax, alpha in zip(axs[::every, 0], alphas[::every]):
    ax.set_ylabel(f'{alpha:1.1f}', rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

for given in givens:

    if movie:
        text = fig.text(0.90, 0.90, f'given: {given:4.2f}', fontsize=fstick, color='grey', ha='right')

    MRT = MRT0*(1.0 - given)

    for rowax, alpha in zip(axs, alphas):
        base = MRT*alpha/(1.0 - alpha)
        for ax, rho in zip(rowax, rhos):
            for line in ax.get_lines():
                line.remove()
            Q = mymodule.Rq*pow(base, 1.0/(rho - 1.0))
            a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
            weq = mymodule.fitness(np.array([a2eq]), np.array([a2eq]), np.array([given]), np.array([alpha]), np.array([rho]))
            w = mymodule.fitness(np.full(numa2, a2eq), a2, np.full(numa2, given), np.full(numa2, alpha), np.full(numa2, rho))
            #a2m = (mymodule.a2max - a2*given*Q*mymodule.b)/(1.0 + Q*mymodule.b*(1.0 - given)) # Optimal a2 given partner's a2
            #a2m[a2m < 0.0] = 0.0
            #a2m[a2m > mymodule.a2max] = mymodule.a2max
            ax.plot(a2, w, linewidth=4, c=cm.magma(weq/traitvmax))

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
