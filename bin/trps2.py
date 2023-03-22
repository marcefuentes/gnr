#! /usr/bin/env python

import os
import time

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

# Options

givens = [1.0, 0.95, 0.5]

num = 21     # Number of subplot rows and columns
nframes = 21 # Number of frames

movie = True
plotsize = 8

# Add data to figure

def update(low, lines):

    lows = np.full(AA.shape, low)
    highs = lows + 0.01
    for g, given in enumerate(givens):

        T = my.fitness(highs, lows, given, AA, RR)
        R = my.fitness(highs, highs, given, AA, RR)
        P = my.fitness(lows, lows, given, AA, RR)
        S = my.fitness(lows, highs, given, AA, RR)
        Z = my.gamecolors(T, R, P, S)
        Ma = np.maximum.reduce([T, R, P, S])
        Mi = np.minimum.reduce([T, R, P, S])
        Tn = (T - Mi)/(Ma - Mi)
        Rn = (R - Mi)/(Ma - Mi)
        Pn = (P - Mi)/(Ma - Mi)
        Sn = (S - Mi)/(Ma - Mi)
        y = np.stack((Tn, Rn, Pn, Sn), axis=-1)
        linecolor = np.full(highs.shape, 'white')
        red = np.full(highs.shape, 'red')
        m = lows > highs
        linecolor[m] = red[m]

        Zg = my.gamecolors(T, R, P, S)
        greys = np.full(Zg.shape, [0.8, 0.8, 0.8, 1.0])
        m = (Zg == my.colormap['white'])
        Zg[m] = greys[m]
        for (a, r, i), _ in np.ndenumerate(y):
            lines[g, a, r].set_ydata(y[a, r])
            lines[g, a, r].axes.set_facecolor(Zg[a, r])
            lcolor = linecolor[a, r] 
            lines[g, a, r].set_color(lcolor)
            lines[g, a, r].set_markerfacecolor(lcolor)

    return lines.flatten()

# Get data

alphas = np.linspace(my.alphamax, my.alphamin, num=num)
logess = np.linspace(my.logesmin, my.logesmax, num=num)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
lows = np.linspace(0.001, 0.999, num=nframes)

# Figure properties

width = plotsize
height = plotsize*len(givens)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
xlim=[0, 5]
ylim=[-0.1, 1.1]
step = int(num/2)
xaxis = [1, 2, 3, 4]
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig = plt.figure(figsize=(width*1.13, height))
outergrid = fig.add_gridspec(nrows=len(givens),
                             ncols=1,
                             left=0.25,
                             right=0.85)
axs = np.empty((len(givens),
                len(alphas),
                len(rhos)),
               dtype=object)

for g, given in enumerate(givens):
    grid = outergrid[g].subgridspec(nrows=num,
                                    ncols=num,
                                    wspace=0,
                                    hspace=0)
    axs[g] = grid.subplots()

left_x = axs[0, 0, 0].get_position().x0
right_x = axs[-1, -1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0, 0].get_position().y1
bottom_y = axs[-1, -1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.2,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=left_x*0.2,
              y=center_y,
              fontsize=biglabels)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    ax.set(xlim=xlim, ylim=ylim)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.1)

letter = ord('a')
for g, given in enumerate(givens):
    axs[g, 0, 0].set_title(chr(letter),
                           fontsize=plotsize*5,
                           pad = 11,
                           weight='bold',
                           loc='left')
    letter += 1
    if g == 0:
        axs[g, 0, 10].set_title('Games',
                                pad=plotsize*9,
                                fontsize=plotsize*5)
    for a in range(0, num, step):
        axs[g, a, 0].set(yticks=[ylim[1]/2])
        axs[g, a, 0].set_yticklabels([alphas[a]],
                                         fontsize=ticklabels)
        for l in range(0, num, step):
            axs[g, -1, l].set(xticks=[xlim[1]/2], xticklabels=[])
            if given == givens[-1]:
                axs[g, -1, l].set_xticklabels([f'{logess[l]:.0f}'],
                                                 fontsize=ticklabels)

# Assign axs objects to variables
# (Line2D objects to lines)

lines = np.empty(axs.shape, dtype=object)
dummy_y = np.zeros_like(xaxis)

for g, given in enumerate(givens):
    for a, alpha in enumerate(alphas):
        for r, rho in enumerate(rhos):
            ax = axs[g, a, r]
            lines[g, a, r], = ax.plot(xaxis,
                                      dummy_y,
                                      linewidth=3,
                                      marker='o',
                                      markerfacecolor='white',
                                      markersize=plotsize/3)

# Add data and save figure

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=lows,
                        fargs=(lines,),
                        blit=True)
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    update(lows[-1], lines,)
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
