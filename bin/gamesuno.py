#! /usr/bin/env python

import matplotlib.pyplot as plt
import mymodule
import numpy as np
import time

start_time = time.perf_counter ()

given = 0.95
alpha = 0.9
loges = 3.0
numa2 = 1000
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
mina = 0.0
X, Y = np.meshgrid(np.linspace(mina, mymodule.a2max, num=numa2),
                    np.linspace(mymodule.a2max, mina, num=numa2))
H = np.copy(Y)
Y[(X > Y)] = X[(X > Y)] 
X[(X > H)] = H[(X > H)] 
RRR, AAA = np.meshgrid(np.repeat(rho, numa2),
                        np.repeat(alpha, numa2))

xmin = mina
xmax = mymodule.a2max
xlabel = '$\it{a}$'
ymin = mina
ymax = mymodule.a2max
ylabel = '$\it{a}$'

xticklabels = [round(xmin, 2),
                round((xmin + xmax)/2, 2),
                round(xmax, 2)]
yticklabels = [round(ymin, 2),
                round((ymin + ymax)/2, 2),
                round(ymax, 2)]
extentg = 0, numa2, 0, numa2

fig, ax = plt.subplots(figsize=(6, 6))

ax.set(xticks=[0, numa2/2, numa2],
        yticks=[0, numa2/2, numa2])
ax.set_xticklabels(xticklabels, fontsize=fstick)
ax.set_yticklabels(yticklabels, fontsize=fstick) 

#if givens[0] > 0.9999999:
#    given = 0.9999999
#for axrow, given in zip(axs, givens):

T = mymodule.fitness(Y, X, given, AAA, RRR)
S = mymodule.fitness(X, Y, given, AAA, RRR)
R = mymodule.fitness(Y, Y, given, AAA, RRR)
P = mymodule.fitness(X, X, given, AAA, RRR)
Z = np.full([numa2, numa2, 4], mymodule.colormap['white'])
mymodule.gamecolors(T, R, P, S, Z)
ax.imshow(Z, extent=extentg)
ax.scatter(a2eq*numa2, a2eq*numa2, s=50, edgecolor='red', facecolor='none', marker='o')
ax.scatter(a2social*numa2, a2social*numa2, s=50, edgecolor='green', facecolor= 'none', marker='o')

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
