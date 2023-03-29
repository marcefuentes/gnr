#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

# Options

folders = ['given95', 'given00']
movie = True
nframes = 21 # Number of frames
plotsize = 8

# Add data to figure

def update(distance, artists):

    lows = a2eqs[0] + distance*(a2eqs[1] - a2eqs[0])
    highs = a2eqs[0] + (distance + 1.0/nframes)*(a2eqs[1] - a2eqs[0])
    T = my.fitness(highs, lows, given, AA, RR)
    R = my.fitness(highs, highs, given, AA, RR)
    P = my.fitness(lows, lows, given, AA, RR)
    S = my.fitness(lows, highs, given, AA, RR)
    Ma = np.maximum.reduce([T, R, P, S])
    Mi = np.minimum.reduce([T, R, P, S])
    Tn = (T - Mi)/(Ma - Mi)
    Rn = (R - Mi)/(Ma - Mi)
    Pn = (P - Mi)/(Ma - Mi)
    Sn = (S - Mi)/(Ma - Mi)
    y = np.stack((Tn, Rn, Pn, Sn), axis=-1)

    for (a, r, i), _ in np.ndenumerate(y):
        artists[a, r].set_ydata(y[a, r])

    return artists.flatten()

# Data

a2eqs = np.empty(len(folders), dtype=object)
for i, folder in enumerate(folders):
    filelist = glob(os.path.join(folder, 'none', '*.csv'))
    df = [my.read_file(file, False) for file in filelist]
    df = pd.concat(df, ignore_index=True)
    if i == 0:
        given = df.Given.iloc[0]
    a = pd.pivot_table(df,
                       values='a2Seenmean',
                       index='alpha',
                       columns='logES')
    a = a.sort_index(axis=0, ascending=False)
    a2eqs[i] = a.to_numpy()

alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1.0 - 1.0/pow(2.0, logess)
RR, AA = np.meshgrid(rhos, alphas)

distances = np.linspace(0.0, 1.0 - 1.0/nframes, num=nframes)

# Figure properties

width = plotsize
height = plotsize
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3.5
xlim=[0, 5]
ylim=[-0.1, 1.1]
step = int(nc/2)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig = plt.figure(figsize=(width, height))
grid = fig.add_gridspec(nrows=nr,
                        ncols=nc,
                        wspace=0,
                        hspace=0,
                        left=0.2,
                        right=0.9,
                        top=0.9,
                        bottom=0.2)
axs = grid.subplots()

left_x = axs[0, 0].get_position().x0
right_x = axs[-1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0].get_position().y1
bottom_y = axs[-1, -1].get_position().y0
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

for i in range(0, nc, step):
    axs[i, 0].set(yticks=[ylim[1]/2])
    axs[i, 0].set_yticklabels([alphas[i]],
                              fontsize=ticklabels)
    axs[-1, i].set(xticks=[xlim[1]/2])
    axs[-1, i].set_xticklabels([f'{logess[i]:.0f}'],
                               fontsize=ticklabels)

# Assign axs objects to variables
# (Line2D)

artists = np.empty(axs.shape, dtype=object) 
xaxis = [1, 2, 3, 4]
dummy_y = np.zeros_like(xaxis)
frames = distances
frame0 = distances[0]

for a, alpha in enumerate(alphas):
    for r, rho in enumerate(rhos):
        ax = axs[a, r]
        artists[a, r], = ax.plot(xaxis,
                                 dummy_y,
                                 linewidth=1,
                                 marker='o',
                                 markersize=plotsize/3)

# Add data and save figure

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(artists,),
                        blit=True)
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    update(frame0, artists,)
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
