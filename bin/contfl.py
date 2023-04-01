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
folders = ['given100', 'given95', 'given50']
subfolders = ['p', 'r']

numa2 = 64
theory = False
movie = False
plotsize = 8

# Add data to figure

def update(t, lines):

    for f, folder in enumerate(folders):
        given = dfprivates[f].Given.iloc[0]
        if theory:
            a2privates = my.a2eq(given, AA, RR)
            ws = my.fitness(a2privates, a2privates, given, AA, RR)
        else:
            a2privates = my.getZ(t, dfprivates[f], 'a2Seenmean')
            ws = my.getZ(t, dfprivates[f], 'wmean')

        for a, alpha in enumerate(alphas):
            for l, rho in enumerate(rhos):
                w = ws[a, l]
                a2s = np.full(xaxis.shape, a2privates[a, l])
                T = my.fitness(a2s, xaxis, given, alpha, rho)
                R = my.fitness(a2s, a2s, given, alpha, rho)
                P = my.fitness(xaxis, xaxis, given, alpha, rho)
                S = my.fitness(xaxis, a2s, given, alpha, rho)
                Ti = my.fitness(xaxis, a2s, given, alpha, rho)
                Ri = my.fitness(xaxis, xaxis, given, alpha, rho)
                Pi = my.fitness(a2s, a2s, given, alpha, rho)
                Si = my.fitness(a2s, xaxis, given, alpha, rho)

                m = xaxis > a2s
                T[m] = Ti[m]
                R[m] = Ri[m]
                P[m] = Pi[m]
                S[m] = Si[m]

                z = my.gamecolors(T, R, P, S)
                y = np.full(xaxis.shape, w)
                lines[f, 0, a, l].set_offsets(np.column_stack((xaxis, y)))
                lines[f, 0, a, l].set_color(z)
                lines[f, 0, a, l].set_sizes([6])

                y = my.fitness(xaxis, xaxis, given, alpha, rho)
                m = xaxis < a2s
                y[m] = -1.0
                m = y < w
                co = np.full(xaxis.shape, 'white')
                co[m] = orange[m]
                lines[f, 1, a, l].set_offsets(np.column_stack((xaxis, y)))
                lines[f, 1, a, l].set_color(co)

                y = my.fitness(xaxis, xaxis, given, alpha, rho) 
                m = y < w
                co = np.full(xaxis.shape, 'white')
                co[m] = orange[m]
                lines[f, 2, a, l].set_offsets(np.column_stack((xaxis, y)))
                lines[f, 2, a, l].set_color(co)

        for c, trait in enumerate(traits):
            Z = my.getZ(t, dftraits[f, c], trait)
            if 'Grain' in trait:
                Z = 1.0 - Z
            for (a, l), _ in np.ndenumerate(Z):
                bgcolor = cm.viridis(Z[a, l]/my.a2max)
                lines[f, c + 1, a, l].axes.set_facecolor(bgcolor)
    if movie:
        fig.texts[2].set_text(f't\n{t}')

    return lines.flatten()

# Data

filelist = glob(os.path.join('given00', 'none', '*.csv'))
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
xlim = [0.0, my.a2max]
ylim = [0.0, my.wmax]
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
             f't\n0',
             fontsize=biglabels,
             color='grey',
             ha='right')

# Assign axs objects to variables
# (PathCollection)

lines = np.empty_like(axs) 
xaxis = np.linspace(0.0, my.a2max, num=numa2)
orange = np.full(xaxis.shape, 'red')
dummy_y = np.zeros_like(xaxis)
frames = ts
frame0 = ts[-1]

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        for a, alpha in enumerate(alphas):
            for l, loges in enumerate(logess):
                ax = axs[f, c, a, l] 
                lines[f, c, a, l] = ax.scatter(xaxis,
                                                 dummy_y,
                                                 marker='s',
                                                 color='white',
                                                 s=1)

# Add data and save figure

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(lines,),
                        blit=True)
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    update(frame0, lines,)
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
