#! /usr/bin/env python

from matplotlib import cm
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import time

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

titles = ['Games',
          '$\it{R}$ - $\it{P}$',
          '$\it{T}$ + $\it{S}$ - 2$\it{R}$']
given = 0.95
num = 5    # Number of subplot rows & columns
distances = [0.2, 0.5, 0.8]

plotsize = 6

alphas = np.linspace(mymodule.alphamax, mymodule.alphamin, num=num)
logess = np.linspace(mymodule.logesmin, mymodule.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)

xlim=[0.0, mymodule.a1max]
ylim=[0.0, mymodule.a2max]
step = int(num/2)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
letterposition = 1.035
xmin = logess[0]
xmax = logess[-1]
ymin = alphas[-1]
ymax = alphas[0]
markersize = 20.0
width = plotsize*len(titles)
height = plotsize
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig = plt.figure(figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.502,
              y=0.0,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=0.07,
              y=0.5,
              fontsize=biglabels)

outergrid = fig.add_gridspec(nrows=1,
                             ncols=len(titles),
                             left=0.15,
                             right=0.85,
                             top=0.8,
                             bottom=0.2)

axss = []
for g, title in enumerate(titles):
    grid = outergrid[g].subgridspec(nrows=num,
                                    ncols=num,
                                    wspace=0,
                                    hspace=0)
    axs = grid.subplots()
    axs[0, int(num/2)].set_title(title,
                                 pad=plotsize*5,
                                 fontsize=plotsize*5)
    axs[0, 0].set_title(chr(letter),
                        fontsize=plotsize*5,
                        weight='bold',
                        loc='left')
    letter += 1

    for ax in fig.get_axes():
        ax.set(xticks=[], yticks=[])
        ax.set(xlim=xlim, ylim=ylim)
    if g == 0:
        for i in range(0, num, step):
            axs[i, 0].set_ylabel(f'{alphas[i]:3.1f}',
                                 rotation='horizontal',
                                 horizontalalignment='right',
                                 verticalalignment='center',
                                 fontsize=ticklabels)
    for j in range(0, num, step):
        axs[-1, j].set_xlabel(f'{logess[j]:2.0f}',
                              fontsize=ticklabels)

    axss.append(axs)

for i, alpha in enumerate(alphas):
    AA = np.array([alpha, alpha, alpha])
    for j, rho in enumerate(rhos):
        RR = np.array([rho, rho, rho])
        a2eq = mymodule.a2eq(given, AA, RR)
        a2social = mymodule.a2eq(0.0, AA, RR)
        X = distances*a2eq
        Y = a2social + distances*(mymodule.a2max - a2social)
        
        T = mymodule.fitness(Y, X, given, AA, RR)
        R = mymodule.fitness(Y, Y, given, AA, RR)
        P = mymodule.fitness(X, X, given, AA, RR)
        S = mymodule.fitness(X, Y, given, AA, RR)

        Z = mymodule.gamecolors(T, R, P, S)
        axss[0][i][j].scatter(X, Y,
                            marker='o',
                            s=markersize,
                            color=Z)
        Z = mymodule.nodilemmacolorsg(T, R, P, S)
        axss[0][i][j].scatter(X, Y,
                            marker='o',
                            s=markersize,
                            color=Z)

        N = mymodule.nodilemmacolors(T, R, P, S)

        Z = R - P
        axss[1][i][j].scatter(X, Y,
                           marker='o',
                           s=markersize,
                           c=cm.viridis(Z))
        axss[1][i][j].scatter(X, Y,
                           marker='o',
                           s=markersize*0.7,
                           c=N)

        Z = T + S - 2.0*R
        axss[2][i][j].scatter(X, Y,
                           marker='o',
                           s=markersize,
                           c=cm.viridis(Z))
        axss[2][i][j].scatter(X, Y,
                           marker='o',
                           s=markersize*0.7,
                           c=N)

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
