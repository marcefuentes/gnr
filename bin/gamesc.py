#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule
import numpy as np
import time

start_time = time.perf_counter ()

traitlabels = ['Games',
                '$\it{T}$ + $\it{S}$ - 2$\it{R}$',
                '$\it{S}$ - $\it{P}$']
given = 0.95
filename = 'gamesc'

num = 1001
numg = 21
numa2 = 64

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
#a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)
MRT = MRT0*(1.0 - given)
Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
weq = mymodule.fitness(a2eq, a2eq, given, AA, RR)

X, Y = np.meshgrid(np.linspace(0.0, mymodule.a2max, num=numa2),
                    np.linspace(mymodule.a2max, 0.0, num=numa2))
X = np.tile(A=X, reps=[numg, numg])
Y = np.tile(A=Y, reps=[numg, numg])
H = np.copy(Y)
Y[(X > Y)] = X[(X > Y)] 
X[(X > H)] = H[(X > H)] 
alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=numg)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=numg)
rhos = 1.0 - 1.0/pow(2, logess)
RRR, AAA = np.meshgrid(np.repeat(rhos, numa2),
                        np.repeat(alphas, numa2))
T = mymodule.fitness(Y, X, given, AAA, RRR)
R = mymodule.fitness(Y, Y, given, AAA, RRR)
P = mymodule.fitness(X, X, given, AAA, RRR)
S = mymodule.fitness(X, Y, given, AAA, RRR)

xmin = logess[0]
xmax = logess[-1]
xlabel = 'Substitutability of $\it{B}$'
ymin = alphas[-1]
ymax = alphas[0]
ylabel = 'Value of $\it{B}$'

xticklabels = [round(xmin),
                round((xmin + xmax)/2),
                round(xmax)]
yticklabels = [round(ymin, 1),
                round((ymin + ymax)/2, 1),
                round(ymax, 1)]
extent2 = 0, num, 0, num
extentg = 0, numg*numa2, 0, numg*numa2

fig, axs = plt.subplots(nrows=1,
                        ncols=len(traitlabels),
                        figsize=(6*len(traitlabels), 6))
fig.supxlabel(xlabel, x=0.513, y=0.01, fontsize=fslabel*1.2)
fig.supylabel(ylabel, x=0.03, y=0.493, fontsize=fslabel*1.2)

letter = ord('a')
for ax, traitlabel in zip(axs, traitlabels):
    ax.text(0, 
            numg*numa2*1.035,
            chr(letter),
            fontsize=fslabel*0.8,
            weight='bold')
    letter += 1
    ax.set_title(traitlabel, pad=40.0, fontsize=fslabel*0.9)
    ax.set(xticks=[0, numg*numa2/2, numg*numa2],
            yticks=[0, numg*numa2/2, numg*numa2],
            xticklabels=[])
    ax.set_xticklabels(xticklabels, fontsize=fstick)
    if ax.get_subplotspec().is_first_col():
        ax.set_yticklabels(yticklabels, fontsize=fstick) 
    else:
        ax.set_yticklabels([]) 

Z = np.full([numg*numa2, numg*numa2, 4], mymodule.colormap['white'])
mymodule.gamecolors(T, R, P, S, Z)
axs[0].imshow(Z, extent=extentg)

Z = np.zeros([numg*numa2, numg*numa2])
mask = (mymodule.prisoner(T, R, P, S) & (2.0*R < T + S)) | (mymodule.deadlock(T, R, P, S) & (2.0*P < T + S))
Z[mask] = 0.1 - (2.0*R[mask] - T[mask] - S[mask])
Z = np.ma.masked_where(Z == 0.0, Z)
cmap = plt.cm.viridis
cmap.set_bad(color='white')
axs[1].imshow(Z, extent=extentg, cmap=cmap)

Z = np.zeros([numg*numa2, numg*numa2])
mask = mymodule.prisoner(T, R, P, S)
Z[mask] = 1.0 - (S[mask] - P[mask])
Z = np.ma.masked_where(Z == 0.0, Z)
cmap = plt.cm.viridis
cmap.set_bad(color='white')
axs[2].imshow(Z, extent=extentg, cmap=cmap)

#axs[2].imshow(a2eq, extent=extent2, cmap='viridis', vmin=0, vmax=traitvmaxs[0])
#axs[3].imshow(weq, extent=extent2, cmap='viridis', vmin=0, vmax=traitvmaxs[1])

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
