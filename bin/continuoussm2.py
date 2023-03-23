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
titles = ['Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner']
folders = ['given100', 'given95', 'given50']
subfolders = ['none', 'p', 'r']

movie = False
plotsize = 8

# Add data to figure

def init(lines):

    highs = pd.pivot_table(dfsocial,
                          values='a2Seenmean',
                          index='alpha',
                          columns='logES')
    highs = highs.sort_index(axis=0, ascending=False)
    highs = highs.to_numpy()

    for f, folder in enumerate(folders):
        df = dfs[f, 0]
        if movie:
            m = df.Time == t
            df = df.loc[m]
        given = df.Given.iloc[0]
        lows = pd.pivot_table(df,
                              values='a2Seenmean',
                              index='alpha',
                              columns='logES')
        lows = lows.sort_index(axis=0, ascending=False)
        lows = lows.to_numpy()
        T = my.fitness(highs, lows, given, AA, RR)
        R = my.fitness(highs, highs, given, AA, RR)
        P = my.fitness(lows, lows, given, AA, RR)
        S = my.fitness(lows, highs, given, AA, RR)

        y = 2.0*R - P - T - S
        ymax = np.max(y)
        ymin = np.min(y)
        y0 = (y - ymin)/(ymax - ymin)

        y = T + S - R - P
        #m = R < P
        #y[m] = T[m] + S[m] - 2.0*P[m]
        ymax = np.max(y)
        ymin = np.min(y)
        y1 = (y - ymin)/(ymax - ymin)

        m = ((R > P) & (P < S)) | ((R < P) & (R < T)) 
        transp = np.ones_like(y)
        transp0 = np.zeros_like(y)
        transp[m] = transp0[m]

        for a, alpha in enumerate(alphas):
            for r, rho in enumerate(rhos):
                lines[f, 0, a, r].set_ydata([y0[a, r], y0[a, r]])
                lines[f, 1, a, r].set_ydata([y1[a, r], y1[a, r]])
                lines[f, 1, a, r].set_alpha(transp[a, r])

    return lines.flatten()

def update(t, lines):
        
    for f, folder in enumerate(folders):
        for c, trait in enumerate(traits):
            df = dfs[f, c + 1]
            if movie:
                m = df.Time == t
                df = df.loc[m]
            Z = pd.pivot_table(df,
                               values=trait,
                               index='alpha',
                               columns='logES')
            Z = Z.sort_index(axis=0, ascending=False)
            Z = Z.to_numpy()
            Z = 1.0 - Z

            for (a, r), _ in np.ndenumerate(Z):
                bgcolor = cm.viridis(Z[a, r]/my.a2max)
                lines[f, c, a, r].axes.set_facecolor(bgcolor)

    return lines.flatten()

# Get data

def read_file(file, alltimes):
    df = pd.read_csv(file)
    if not alltimes:
        df = df.tail(1)
    return df

filelist = glob('given00/none/*.csv')
df = [read_file(file, False) for file in filelist]
dfsocial = pd.concat(df, ignore_index=True)

dfs = np.empty((len(folders), len(subfolders)), dtype=object)
for i, folder in enumerate(folders):
    for j, subfolder in enumerate(subfolders):
        filelist = glob(os.path.join(folder, subfolder, '*.csv'))
        df = [read_file(file, movie) for file in filelist]
        dfs[i, j] = pd.concat(df, ignore_index=True)

df = dfs[0, 0]
ts = df.Time.unique()
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1.0 - 1.0/pow(2.0, logess)
RR, AA = np.meshgrid(rhos, alphas)

# Figure properties

width = plotsize*len(traits)
height = plotsize*len(folders)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
xlim = [0, 1]
ylim=[-0.1, 1.1]
step = int(nr/2)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=len(folders),
                             ncols=len(traits))
axs = np.empty((len(folders),
                len(traits),
                nr,
                nc),
                dtype=object)

for f, folder in enumerate(folders):
    for c, trait in enumerate(traits):
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
for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        axs[f, c, 0, 0].set_title(chr(letter),
                                  fontsize=plotsize*5,
                                  pad = 11,
                                  weight='bold',
                                  loc='left')
        letter += 1
        if f == 0:
            axs[f, c, 0, 10].set_title(title,
                                       pad=plotsize*9,
                                       fontsize=plotsize*5)
        for a in range(0, nr, step):
            axs[f, c, a, 0].set(yticks=[ylim[1]/2], yticklabels=[])
            if c == 0:
                axs[f, c, a, 0].set_yticklabels([alphas[a]],
                                                fontsize=ticklabels)
        for l in range(0, nc, step):
            axs[f, c, -1, l].set(xticks=[xlim[1]/2], xticklabels=[])
            if folder == folders[-1]:
                axs[f, c, -1, l].set_xticklabels([f'{logess[l]:.0f}'],
                                                 fontsize=ticklabels)

# Assign axs objects to variables
# (Line2D objects to lines)

lines = np.empty(axs.shape, dtype=object)
dummy_y = np.zeros_like(xlim)

for f, folder in enumerate(folders):
    for c, trait in enumerate(traits):
        for a, alpha in enumerate(alphas):
            for r, rho in enumerate(rhos):
                ax = axs[f, c, a, r] 
                lines[f, c, a, r], = ax.plot(xlim,
                                             dummy_y,
                                             linewidth=1,
                                             color='white')

# Add data and save figure

init(lines,)

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=ts,
                        fargs=(lines,),
                        blit=True)
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    update(ts[-1], lines,)
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
