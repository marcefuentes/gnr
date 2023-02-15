#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule
import numpy as np
import time

start_time = time.perf_counter ()

traitlabels = ['Games', 'Effort to get $\it{B}$', 'Fitness']
given = 0.95

numg = 21
numa2 = 64
filename = 'gamesc'

fslabel = 32 # Label font size
fstick = 18 # Tick font size
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

if given > 0.9999999:
    given = 0.9999999
alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=numg)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=numg)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
MRT0 = mymodule.b*mymodule.Rq
Q0 = mymodule.Rq*pow(MRT0*AA/(1.0 - AA), 1.0/(RR - 1.0))
#a2social = mymodule.a2max/(1.0 + Q0*mymodule.b)
MRT = MRT0*(1.0 - given)
Q = mymodule.Rq*pow(MRT*AA/(1.0 - AA), 1.0/(RR - 1.0))
a2eq = mymodule.a2max/(1.0 + Q*mymodule.b)
weq = mymodule.fitness(a2eq, a2eq, given, AA, RR)
a2eq = a2eq*0.0
X = np.zeros([numg*numa2, numg*numa2])
Y = np.zeros([numg*numa2, numg*numa2])
for i in range(0, numg):
    for j in range(0, numg):
        for k in range(numa2):
            for l in reversed(range(numa2)):
                X[i*numa2+k, j*numa2+l] = a2eq[i, j] + (1.0 - a2eq[i, j])*k/numa2
                Y[i*numa2+k, j*numa2+l] = a2eq[i, j] + (1.0 - a2eq[i, j])*l/numa2
H = np.copy(Y)
Y[(X > Y)] = X[(X > Y)] 
X[(X > H)] = H[(X > H)] 
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

xticklabels = [round(xmin),
                round((xmin + xmax)/2),
                round(xmax)]
yticklabels = [round(ymin, 1),
                round((ymin + ymax)/2, 1),
                round(ymax, 1)]
extentg = 0, numg*numa2, 0, numg*numa2

fig, ax = plt.subplots(figsize=(6, 6))
fig.supxlabel(xlabel, x=0.513, y=0.01, fontsize=fslabel*1.2)
fig.supylabel(ylabel, x=0.03, y=0.493, fontsize=fslabel*1.2)

ax.set(xticks=[0, numg*numa2/2, numg*numa2],
        yticks=[0, numg*numa2/2, numg*numa2])
ax.set_yticklabels(yticklabels, fontsize=fstick) 
ax.set_xticklabels(xticklabels, fontsize=fstick)

ax.imshow(Z, extent=extentg)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
