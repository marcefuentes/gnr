#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule
import numpy as np
import time

start_time = time.perf_counter ()

traitlabels = ['Games',
                '$\it{T}$ + $\it{S}$ - 2$\it{R}$',
                '$\it{T}$ + $\it{S}$']
given = 0.95
alpha = 0.6
loges = 0.5
numa2 = 1025
filename = 'gamesuno'

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

rho = 1.0 - 1.0/pow(2, loges)
MRT0 = mymodule.b*mymodule.Rq
Q0 = mymodule.Rq*pow(MRT0*alpha/(1.0 - alpha), 1.0/(rho - 1.0))
a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)
MRT = MRT0*(1.0 - given)
Q = mymodule.Rq*pow(MRT*alpha/(1.0 - alpha), 1.0/(rho - 1.0))
a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)

X, Y = np.meshgrid(np.linspace(0.0, mymodule.a2max, num=numa2),
                    np.linspace(mymodule.a2max, 0.0, num=numa2))
H = np.copy(Y)
Y[(X > Y)] = X[(X > Y)] 
X[(X > H)] = H[(X > H)] 
RRR, AAA = np.meshgrid(np.repeat(rho, numa2),
                        np.repeat(alpha, numa2))
T = mymodule.fitness(Y, X, given, AAA, RRR)
R = mymodule.fitness(Y, Y, given, AAA, RRR)
P = mymodule.fitness(X, X, given, AAA, RRR)
S = mymodule.fitness(X, Y, given, AAA, RRR)

xmin = 0.0
xmax = mymodule.a2max
xlabel = 'Effort to get $\it{B}$'
ymin = 0.0
ymax = mymodule.a2max
ylabel = 'Effort to get $\it{B}$'

xticklabels = [round(xmin, 2),
                round((xmin + xmax)/2, 2),
                round(xmax, 2)]
yticklabels = [round(ymin, 2),
                round((ymin + ymax)/2, 2),
                round(ymax, 2)]
extentg = 0, numa2, 0, numa2

fig, axs = plt.subplots(nrows=1,
                        ncols=len(traitlabels),
                        figsize=(6*len(traitlabels), 6))
fig.supxlabel(xlabel, x=0.513, y=0.01, fontsize=fslabel*1.2)
fig.supylabel(ylabel, x=0.03, y=0.493, fontsize=fslabel*1.2)

letter = ord('a')
for ax, traitlabel in zip(axs, traitlabels):
    ax.set(yticks=[0, numa2/2, numa2])
    ax.set(xticks=[0, numa2/2, numa2])
    ax.set_title(traitlabel, pad=40.0, fontsize=fslabel*0.9)
    ax.set_xticklabels(xticklabels, fontsize=fstick)
    ax.text(0, 
            numa2*1.035,
            chr(letter),
            fontsize=fslabel*0.8,
            weight='bold')
    letter += 1
    if ax.get_subplotspec().is_first_col():
        ax.set_yticklabels(yticklabels, fontsize=fstick) 
    else:
        ax.set_yticklabels([])

if given > 0.9999999:
    given = 0.9999999

Z = np.full([numa2, numa2, 4], mymodule.colormap['white'])
mymodule.gamecolors(T, R, P, S, Z)
axs[0].imshow(Z, extent=extentg)
axs[0].scatter(a2eq*numa2, a2eq*numa2, s=50, edgecolor='red', facecolor='none', marker='o')
axs[0].scatter(a2social*numa2, a2social*numa2, s=50, edgecolor='green', facecolor= 'none', marker='o')

Z = np.zeros([numa2, numa2])
mask = (mymodule.prisoner(T, R, P, S) & (2.0*R < T + S)) | (mymodule.deadlock(T, R, P, S) & (2.0*P < T + S))
Z[mask] = 0.1 - (2.0*R[mask] - T[mask] - S[mask])
Z = np.ma.masked_where(Z == 0.0, Z)
cmap = plt.cm.viridis
cmap.set_bad(color='white')
axs[1].imshow(Z, extent=extentg, cmap=cmap)

Z = np.zeros([numa2, numa2])
Z = T + S
vmax = mymodule.fitness(np.array([mymodule.a2max]),
                                np.array([mymodule.a2max]),
                                np.array([0.0]),
                                np.array([0.9]),
                                np.array([5.0]))
axs[2].imshow(Z, extent=extentg, cmap=cmap, vmin=0, vmax=vmax)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
