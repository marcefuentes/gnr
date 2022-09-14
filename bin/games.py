#! /usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import time

start_time = time.perf_counter ()

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

alpha = 0.5
R1 = 2.0
R2 = 2.0
R = R2/R1
a1max = 1.0
a2max = 1.0
b = a2max/a1max

npoints = 32

num = 11
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.0
maxgiven = 1.0
aC = 0.5
aD = 0.0

log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
givens = np.linspace(maxgiven, mingiven, num=num)
givens[0] = 0.99999
Xrhos, Ygivens = np.meshgrid(rhos, givens)

x = np.linspace(0.001, 0.999, npoints)
a2partner = x
y = np.linspace(0.999, 0.001, npoints)
X, Y = np.meshgrid(x, y)

def a2eq(given, rho):
    T = b*R*(1.0 - given)
    Q = R*pow(T*(1.0 - alpha)/alpha, 1.0/(rho - 1.0))
    a2 = a2max/(1.0 + Q*b)
    return a2

def fitness(x, y, given, rho):
    q1 = (1.0 - y)*R1
    q2 = y*R2*(1.0 - given) + x*R2*given
    if rho == 0.0:
        w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

def fitness2(x, y, given, rho):
    q1 = (1.0 - x)*R1
    q2 = x*R2*(1.0 - given) + y*R2*given
    if rho == 0.0:
        w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

def a2maxw(x, given, rho):
    T = b*R*(1.0 - given)
    Q = R*pow(T*(1.0 - alpha)/alpha, 1.0/(rho - 1.0))
    a2 = (a2max - x*given*Q*b)/(1.0 + Q*b*(1.0 - given))
    a2 = 1.0 - a2
    return a2

fslabel = 26 # Label font size
fstick = 18 # Tick font size

fig = plt.figure() 

gs0 = gridspec.GridSpec(1, 2, figure=fig)

gs00 = gridspec.GridSpecFromSubplotSpec(num, num, subplot_spec=gs0[0])
axs = gs00.subplots()

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        Z = fitness(X, Y, given, rho)
        ax.imshow(Z, cmap='magma', vmin=0, vmax=2)
        xaxis = a2partner*npoints
        yaxis = a2maxw(a2partner, given, rho)*npoints
        #ax.plot(npoints - yaxis, npoints - xaxis, color='grey')
        ax.plot(xaxis, yaxis, color='white', alpha=0.6)
        #ypoint = a2eq(given, rho)*npoints
        #ax.scatter(ypoint, npoints - ypoint, color='white', s=0.3)
        ax.set(xticks=[], yticks=[], xlim=(0, npoints-1), ylim=(npoints-1, 0))
        if (given == 0.5) and (rho == rhos[0]):
            ax.set_ylabel("Partner's share of $\it{A}$", fontsize=fslabel)
for ax, log_es in zip(axs[-1, ::5], log_ess[::5]):
    ax.set_xlabel(round(log_es), fontsize=fstick)
for ax, given in zip(axs[::5, 0], givens[::5]):
    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)

#subfigs[0].suptitle('a', fontsize=fslabel, horizontalalignment='left', weight='bold')
#subfigs[0].subplots_adjust(wspace=0, hspace=0)

gs01 = gridspec.GridSpecFromSubplotSpec(num, num, subplot_spec=gs0[1])
axs = gs01.subplots()

xaxis = [1, 2, 3, 4]
givens[0] = 1.0

for row, given in zip(axs, givens):
    for ax, rho in zip(row, rhos):
        R = wC1 = fitness2(aC, aC, given, rho)
        P = wD0 = fitness2(aD, aD, given, rho)
        T = wD1 = fitness2(aD, aC, given, rho)
        S = wC0 = fitness2(aC, aD, given, rho)
        yaxis = [T, R, P, S]
        if (T<R) and (P<S):
            rgb = (0.5, 0.0, 1.0)
        else:
            if (T>R) and (P<S):
                rgb = (0.0, 1.0, 0.0) # Snowdrift
            else:
                rgb = (1.0, 0.0, 0.5)
        ax.plot(xaxis, yaxis, color=rgb, marker='o', markerfacecolor='white', linewidth=1.0, markersize=3)
        ax.set(xticks=[], yticks=[], xlim=(0, 5), ylim=(0.0, 2.0))

for ax, log_es in zip(axs[-1, ::5], log_ess[::5]):
    ax.set_xlabel(round(log_es), fontsize=fstick)

#gs[1].subplots_adjust(wspace=0, hspace=0)

plt.suptitle("Substitutability of $\it{A}$", fontsize=fslabel)

plt.savefig('games.png', dpi=100)
plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
