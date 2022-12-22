#! /usr/bin/env python

import os
import imageio.v2 as iio
import matplotlib.pyplot as plt
import numpy as np
import time

start_time = time.perf_counter ()

minalpha = 0.1
maxalpha = 0.9
minlog_es = -5.0
maxlog_es = 5.0
mingiven = 0.95
maxgiven = 0.95

movie = True

num = 21    # Number of subplot rows and columns
numa2 = 64
filename = 'games_c_es_alpha'
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0

traitlabels = ['Effort to get $\it{B}$', 'Fitness']
traitvmaxs = [1.0, 2.0, 1.8]
fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def fitness(x, y, given, alpha, rho):
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - given) + x*R2*given
    w = q1*q2
    mask = (w > 0.0) & (rho == 0.0)
    w[mask] = pow(q1[mask], 1.0 - alpha[mask])*pow(q2[mask], alpha[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = (1.0 - alpha[mask])*pow(q1[mask], rho[mask]) + alpha[mask]*pow(q2[mask], rho[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = pow(w[mask], 1.0/rho[mask])
    mask = (rho > 0.0)
    w[mask] = pow((1.0 - alpha[mask])*pow(q1[mask], rho[mask]) + alpha[mask]*pow(q2[mask], rho[mask]), 1.0/rho[mask])
    return w

if movie:
    givens = np.linspace(0.0, 1.0, num=num)
    frames = []
else:
    givens = np.array([mingiven])
if givens[0] > 0.9999999:
    givens[0] = 0.9999999
log_ess = np.linspace(minlog_es, maxlog_es, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
alphas = np.linspace(maxalpha, minalpha, num=num)
RR, AA = np.meshgrid(rhos, alphas)
RRR, AAA = np.meshgrid(np.repeat(rhos, numa2), np.repeat(alphas, numa2))
minx = minlog_es
maxx = maxlog_es
miny = minalpha
maxy = maxalpha
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
nc = num
nr = num

xticklabels = [round(minx), round((minx + maxx)/2), round(maxx)]
yticklabels = [miny, (miny + maxy)/2, maxy]
extent = 0, nr, 0, nc
extenta2 = 0, numa2*nr, 0, numa2*nc
prisoner = np.full((numa2*nr, numa2*nc, 4), [0.5, 0.0, 0.0, 1.0])
snowdrift = np.full((numa2*nr, numa2*nc, 4), [0.0, 1.0, 1.0, 1.0])
nodilemma = np.full((numa2*nr, numa2*nc, 4), [1.0, 1.0, 1.0, 1.0])
green = np.full((numa2*nr, numa2*nc, 4), [0.0, 1.0, 0.0, 1.0])

b = a2max/a1max
Rq = R2/R1
MRT0 = b*Rq
X, Y = np.meshgrid(np.linspace(0.0, a2max, num=numa2), np.linspace(a2max, 0.0, num=numa2))
X = np.tile(A=X, reps=[num, num])
Y = np.tile(A=Y, reps=[num, num])

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 18))
fig.delaxes(axs[0, 1])
fig.supxlabel(xlabel, x=0.515, y=0.03, fontsize=fslabel)
fig.supylabel(ylabel, x=0.04, y=0.510, fontsize=fslabel)

letter = ord('b')
for axrow in axs:
    for ax in axrow:
        if ax.get_subplotspec().is_first_row():
            ax.set(xticks=[0, numa2*nc/2, numa2*nc], yticks=[0, numa2*nr/2, numa2*nr], xticklabels=[], yticklabels=[])
            ax.set_title('Game types', pad=50.0, fontsize=fslabel)
            ax.text(0, numa2*nr*1.035, 'a', fontsize=fslabel, weight='bold')
            pos = ax.get_position()
            newpos = [pos.x0, pos.y0+0.04, pos.width, pos.height]
            ax.set_position(newpos)
        else:
            ax.set(xticks=[0, nc/2, nc], yticks=[0, nr/2, nr], xticklabels=[], yticklabels=[])
            ax.text(0, nr*1.035, chr(letter), fontsize=fslabel, weight='bold')
            letter += 1
        if ax.get_subplotspec().is_last_row():
            ax.set_xticklabels(xticklabels, fontsize=fstick)
        if ax.get_subplotspec().is_first_col():
            ax.set_yticklabels(yticklabels, fontsize=fstick) 
for ax, traitlabel in zip(axs[1], traitlabels):
    ax.set_title(traitlabel, pad=50.0, fontsize=fslabel)

for given in givens:

    if movie:
        text = fig.text(0.80, 0.80, f'given\n{given:4.2f}', fontsize=fstick+4, color='grey', ha='right')

    Z = np.copy(green)
    T = fitness(Y, X, given, AAA, RRR)
    R = fitness(Y, Y, given, AAA, RRR)
    P = fitness(X, X, given, AAA, RRR)
    S = fitness(X, Y, given, AAA, RRR)
    mask = (R < P)
    H = R[mask]
    R[mask] = P[mask]
    P[mask] = H
    H = T[mask]
    T[mask] = S[mask]
    S[mask] = H
    mask = (T > R) & (P > S)
    Z[mask] = prisoner[mask]
    mask = (T >= R) & (P <= S) & (R != P)
    Z[mask] = snowdrift[mask]
    mask = ((T < R) & (P < S)) | (R == P)
    Z[mask] = nodilemma[mask]
    #Z = np.tril(Z, k=-1)
    #Z = np.ma.masked_where(Z == 0.0, Z)

    MRT = MRT0*(1.0 - given)
    Q0 = Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2social = a2max/(1.0 + Q0*b)
    wsocial = fitness(a2social, a2social, given, AA, RR)
    Q = Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
    a2eq = a2max/(1.0 + Q*b)
    weq = fitness(a2eq, a2eq, given, AA, RR)
    Mss = [[a2social, wsocial], [a2eq, weq]]

    axs[0, 0].imshow(Z, extent=extenta2)

    for row, Ms in zip(axs[1:], Mss):
        for ax, M, traitvmax in zip(row, Ms, traitvmaxs):
            ax.imshow(M, extent=extent, cmap='magma', vmin=0, vmax=traitvmax)

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
