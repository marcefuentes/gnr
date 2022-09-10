#! /usr/bin/env python

import time
import matplotlib.pyplot as plt
import numpy as np

start_time = time.perf_counter ()

alpha = 0.5
R1 = 2.0
R2 = 2.0
R = R2/R1
a1max = 1.0
a2max = 1.0
b = a2max/a1max
npoints = 32

log_ess = np.linspace(-5, 5, num=11)
rhos = 1.0 - 1.0/pow(2, log_ess)
givens = np.linspace(1.0, 0.0, num=11)
givens[0] = 0.99999
Xrhos, Ygivens = np.meshgrid(rhos, givens)

x = np.linspace(0.001, 0.999, npoints)
a2partner = x
y = np.linspace(0.999, 0.001, npoints)
X, Y = np.meshgrid(x, y)

def fitness(x, y, given, rho):
    q1 = (1.0 - y)*R1
    q2 = y*R2*(1.0 - given) + x*R2*given
    if rho == 0.0:
        w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

def a2maxw(given, rho):
    T = b*R*(1.0 - given)
    Q = R*pow(T*(1.0 - alpha)/alpha, 1.0/(rho - 1.0))
    a2 = (a2max - a2partner*given*Q*b)/(1.0 + Q*b*(1.0 - given))
    a2 = 1.0 - a2
    return a2

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fslabel = 26 # Label font size
fstick = 18 # Tick font size

fig = plt.figure(figsize=(6, 6))
fig.supylabel("Partner's share of $\it{A}$", x=0.03, y=0.54, fontsize=fslabel)
fig.supxlabel("Substitutability of $\it{A}$", x=0.553, y=0.05, fontsize=fslabel)

grid = fig.add_gridspec(len(givens), len(log_ess))

axs = grid.subplots()
plt.subplots_adjust(hspace=0, wspace=0, bottom=0.2, left=0.2)

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        Z = fitness(X, Y, given, rho)
        ax.imshow(Z, cmap='magma', vmin=0, vmax=2)
        xaxis = a2partner*npoints
        yaxis = a2maxw(given, rho)*npoints
        ax.plot(xaxis, yaxis, color='white')
        ax.set(xticks=[], yticks=[], xlim=(0, npoints-1), ylim=(npoints-1, 0))
for ax, log_es in zip(axs[-1, ::5], log_ess[::5]):
    ax.set_xlabel(round(log_es), fontsize=fstick)
for ax, given in zip(axs[::5, 0], givens[::5]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

plt.savefig('games.png', dpi=100)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
