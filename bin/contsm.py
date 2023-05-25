#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib import cm
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

# Options

traits = ['ChooseGrainmean',
          'MimicGrainmean']
titles = ['Games',
          'Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner']
folders = ['given100', 'given095', 'given050']
subfolders = ['p', 'r']

movie = False
plotsize = 6

# Add data to figure

def init(artists):

    for f, folder in enumerate(folders):
        lows = my.getZ(ts[-1], dfprivates[f], 'a2Seenmean')
        highs = lows + 0.01
        given = dftraits[f, 0].Given.iloc[0]
        T = my.fitness(highs, lows, given, AA, RR)
        R = my.fitness(highs, highs, given, AA, RR)
        P = my.fitness(lows, lows, given, AA, RR)
        S = my.fitness(lows, highs, given, AA, RR)

        CG = (S - P + 0.02)*plotsize*400
        zeros = np.zeros_like(CG)
        #m = my.harmony(T, R, P, S) | my.deadlock(T, R, P, S) | my.snowdrift(T, R, P, S)
        #CG[m] = zeros[m]
        #m = CG < 0.0
        #CG[m] = zeros[m]
        #MG = (P - S)*plotsize*10
        MG = (highs - lows)*plotsize*10 
        m = MG < 0.0
        MG[m] = zeros[m]
        #m = ((R > P) & (P < S)) | ((R < P) & (R < T)) 
        #MG[m] = zeros[m]

        Zg = my.gamecolors(T, R, P, S)
        for a, alpha in enumerate(alphas):
            for l, loges in enumerate(logess):
                artists[f, 0, a, l].axes.set_facecolor(Zg[a, l])
                artists[f, 1, a, l].set_sizes([CG[a, l]])
                artists[f, 2, a, l].set_sizes([MG[a, l]])

    return artists.flatten()

def update(t, artists):
    for f, folder in enumerate(folders):
        for c, trait in enumerate(traits):
            Z = my.getZ(t, dftraits[f, c], trait)
            if 'Grain' in trait:
                Z = 1.0 - Z
            for (a, l), _ in np.ndenumerate(Z):
                bgcolor = cm.viridis(Z[a, l]/my.a2max)
                artists[f, c + 1, a, l].axes.set_facecolor(bgcolor)
    if movie:
        fig.texts[2].set_text(f't\n{t}')
    return artists.flatten()

# Data

filelist = glob(os.path.join('given000', 'none', '*.csv'))
dfsocial = my.read_files(filelist, movie)

dfprivates = np.empty(len(folders), dtype=object)
for f, folder in enumerate(folders):
    filelist = glob(os.path.join(folder, 'none', '*.csv'))
    dfprivates[f] = my.read_files(filelist, movie)

dftraits = np.empty((len(folders), len(subfolders)), dtype=object)
for f, folder in enumerate(folders):
    for c, subfolder in enumerate(subfolders):
        filelist = glob(os.path.join(folder, subfolder, '*.csv'))
        dftraits[f, c] = my.read_files(filelist, movie)

df = dftraits[0, 0]
ts = df.Time.unique()
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1.0 - 1.0/pow(2.0, logess)
RR, AA = np.meshgrid(rhos, alphas)

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(folders)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
xlim = [0, 1]
ylim = [0, 1]
step = int(nr/2)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=len(folders),
                             ncols=len(titles))
axs = np.empty((len(folders),
                len(titles),
                nr,
                nc),
               dtype=object)

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        grid = outergrid[f, c].subgridspec(nrows=nr,
                                           ncols=nc,
                                           wspace=0,
                                           hspace=0)
        axs[f, c] = grid.subplots()

left_x = axs[0, 0, 0, 0].get_position().x0
right_x = axs[-1, -1, -1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0, 0, 0].get_position().y1
bottom_y = axs[-1, -1, -1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.3,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=left_x*0.4,
              y=center_y,
              fontsize=biglabels)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    ax.set(xlim=xlim, ylim=ylim)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.1)

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        letter = ord('a') + f*len(titles) + c
        axs[f, c, 0, 0].set_title(chr(letter),
                                  fontsize=plotsize*5,
                                  pad = 11,
                                  weight='bold',
                                  loc='left')
        if f == 0:
            axs[0, c, 0, 10].set_title(title,
                                       pad=plotsize*9,
                                       fontsize=plotsize*5)
        for a in range(0, nr, step):
            axs[f, c, a, 0].set(yticks=[ylim[1]/2.0], yticklabels=[])
            if c == 0:
                axs[f, 0, a, 0].set_yticklabels([alphas[a]],
                                                fontsize=ticklabels)
        for l in range(0, nc, step):
            axs[f, c, -1, l].set(xticks=[xlim[1]/2.0], xticklabels=[])
            if folder == folders[-1]:
                axs[-1, c, -1, l].set_xticklabels([f'{logess[l]:.0f}'],
                                                 fontsize=ticklabels)
if movie:
    fig.text(right_x,
             bottom_y*0.5,
             't\n0',
             fontsize=biglabels,
             color='grey',
             ha='right')

# Assign axs objects to variables
# (PathCollection)

artists = np.empty_like(axs) 
x = [0.5]
y = [0.5]
dummy_z = [0.0]
frames = ts
frame0 = ts[-1]

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        for a, alpha in enumerate(alphas):
            for l, loges in enumerate(logess):
                ax = axs[f, c, a, l] 
                artists[f, c, a, l] = ax.scatter(x,
                                                 y,
                                                 color='white',
                                                 s=dummy_z)

# Add data and save figure

init(artists,)

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
