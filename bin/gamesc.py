#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule
import numpy as np
import time

start_time = time.perf_counter ()

traitlabels = ['Games', 'Effort to get $\it{B}$', 'Fitness']
given = 0.95

num = 1001
numg = 21
numa2 = 64
filename = 'gamesc'

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
S = mymodule.fitness(X, Y, given, AAA, RRR)
R = mymodule.fitness(Y, Y, given, AAA, RRR)
P = mymodule.fitness(X, X, given, AAA, RRR)
Z = np.full([numg*numa2, numg*numa2, 4], mymodule.colormap['white'])
mymodule.gamecolors(T, R, P, S, Z)

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
extentg = 0, numg*numa2, 0, numg*numa2

fig, axs = plt.subplots(nrows=1,
                        ncols=len(traitlabels),
                        figsize=(6*len(traitlabels), 6))
fig.supxlabel(xlabel, x=0.513, y=0.01, fontsize=fslabel*1.2)
fig.supylabel(ylabel, x=0.03, y=0.493, fontsize=fslabel*1.2)

letter = ord('a')
for ax, traitlabel in zip(axs, traitlabels):
    if ax.get_subplotspec().is_first_col():
        ax.text(0, 
                numg*numa2*1.035,
                chr(letter),
                fontsize=fslabel*0.8,
                weight='bold')
        ax.set(xticks=[0, numg*numa2/2, numg*numa2],
                yticks=[0, numg*numa2/2, numg*numa2],
                xticklabels=[])
        ax.set_yticklabels(yticklabels, fontsize=fstick) 
    else:
        ax.text(0, 
                num*1.035,
                chr(letter),
                fontsize=fslabel*0.8,
                weight='bold')
        ax.set(xticks=[0, num/2, num],
                        yticks=[0, num/2, num],
                        xticklabels=[],
                        yticklabels=[])
    letter += 1
    ax.set_title(traitlabel, pad=40.0, fontsize=fslabel*0.9)
    ax.set_xticklabels(xticklabels, fontsize=fstick)

axs[0].imshow(Z, extent=extentg)
axs[1].imshow(a2eq, extent=extent, cmap='viridis', vmin=0, vmax=traitvmaxs[0])
axs[2].imshow(weq, extent=extent, cmap='viridis', vmin=0, vmax=traitvmaxs[1])

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
