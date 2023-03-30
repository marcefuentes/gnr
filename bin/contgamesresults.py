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

traits = ['ChooseGrainmean',
          'MimicGrainmean']
titles = ['Games',
          'Severity of\nsocial dilemma', 
          'Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner']
vmaxs = [2.0, my.a2max, my.a2max]
folders = ['given100', 'given95', 'given50']
subfolders = ['p', 'r']

movie = False
plotsize = 6
ext = 256

# Add data to figure

def init(artists):

    for f, folder in enumerate(folders):
        given = dfprivates[f].Given.iloc[0]
        for a, alpha in enumerate(alphas):
            AA = np.full(X.shape, alpha)
            for l, (rho, loges) in enumerate(zip(rhos, logess)):
                RR = np.full(X.shape, rho)
                T = my.fitness(Y, X, given, AA, RR)
                R = my.fitness(Y, Y, given, AA, RR)
                P = my.fitness(X, X, given, AA, RR)
                S = my.fitness(X, Y, given, AA, RR)
                Z = my.gamecolors(T, R, P, S)
                Z[X >= Y] = [0.9, 0.9, 0.9, 1.0]
                artists[f, a, l].set_array(Z)
    return artists.flatten()

def update(t, artists):
    Zsocial = getZ(t, dfsocial, 'wmean')
    for f, folder in enumerate(folders):
        Z = getZ(t, dfprivates[f], 'wmean')
        Z = Zsocial - Z       
        artists[f, 0].set_array(Z)
        for c, trait in enumerate(traits):
            Z = getZ(t, dftraits[f, c], trait)
            if 'Grain' in trait:
                Z = 1.0 - Z
            artists[f, c + 1].set_array(Z)
    if movie:
        fig.texts[2].set_text(f't\n{t}')
    return artists.flatten()

def getZ(t, df, trait):
    if movie:
        m = df.Time == t
        df = df.loc[m]
    Z = pd.pivot_table(df,
                       values=trait,
                       index='alpha',
                       columns='logES')
    Z = Z.sort_index(axis=0, ascending=False)
    return Z

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

df = dfsocial
ts = df.Time.unique()
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1.0 - 1.0/pow(2.0, logess)

x = np.linspace(0.0, my.a2max, num=ext)
y = np.flip(x)
X, Y = np.meshgrid(x, y)

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(folders)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
step = int(nr/2)
xticks = [0, nc/2 - 0.5, nc - 1]
yticks = [0, nr/2 - 0.5, nr - 1]
xmin = logess[0]
xmax = logess[-1]
ymin = alphas[-1]
ymax = alphas[0]
xticklabels = [f'{xmin:2.0f}',
               f'{(xmin + xmax)/2.0:2.0f}',
               f'{xmax:2.0f}']
yticklabels = [f'{ymax:3.1f}',
               f'{(ymin + ymax)/2.0:3.1f}',
               f'{ymin:3.1f}']
extent = 0, ext, 0, ext
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=len(folders),
                             ncols=len(titles))
axgames = np.empty((len(folders),
                    nr,
                    nc),
                   dtype=object)
axtraits = np.empty((len(folders),
                     len(titles) - 1), 
                    dtype=object)

for f, folder in enumerate(folders):
    grid = outergrid[f, 0].subgridspec(nrows=nr,
                                       ncols=nc,
                                       wspace=0,
                                       hspace=0)
    axgames[f] = grid.subplots()
    for c in range(len(titles) - 1):
        grid = outergrid[f, c + 1].subgridspec(nrows=1,
                                               ncols=1)
        axtraits[f, c] = grid.subplots()

left_x = axgames[0, 0, 0].get_position().x0
right_x = axtraits[-1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axtraits[0, 0].get_position().y1
bottom_y = axtraits[-1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.3,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=left_x*0.5,
              y=center_y,
              fontsize=biglabels)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])

letter = ord('a')
letterposition = 1.035
for f, folder in enumerate(folders):
    for a, alpha in enumerate(alphas):
        for l, loges in enumerate(logess):
            for axis in ['top','bottom','left','right']:
                axgames[f, a, l].spines[axis].set_linewidth(0.1)
    axgames[f, 0, 0].set_title(chr(letter),
                               pad=plotsize*5/3,
                               fontsize=plotsize*5,
                               weight='bold')
    letter += 1
    if folder == folders[0]:
        axgames[f, 0, 10].set_title(titles[0],
                                    pad=plotsize*9,
                                    fontsize=plotsize*5)
    for a in range(0, nr, step):
        axgames[f, a, 0].set(yticks=[ext/2.0])
        axgames[f, a, 0].set_yticklabels([f'{alphas[a]:.1f}'],
                                         rotation='horizontal',
                                         horizontalalignment='right',
                                         verticalalignment='center',
                                         fontsize=ticklabels)
    for l in range(0, nc, step):
        axgames[f, -1, l].set(xticks=[ext/2.0], xticklabels=[]) 
        if folder == folders[-1]:
            axgames[f, -1, l].set_xticklabels([f'{logess[l]:.0f}'],
                                              fontsize=ticklabels)
    for c in range(len(titles) - 1):
        axtraits[f, c].set(xticks=xticks, yticks=yticks)
        axtraits[f, c].set(xticklabels=[], yticklabels=[])
        axtraits[f, c].text(0,
                            letterposition,
                            chr(letter),
                            transform=axtraits[f, c].transAxes,
                            fontsize=plotsize*5,
                            weight='bold')
        letter += 1
        if folder == folders[0]:
            axtraits[f, c].set_title(titles[c + 1],
                                pad=plotsize*9,
                                fontsize=plotsize*5)
        if folder == folders[-1]:
            axtraits[f, c].set_xticklabels(xticklabels,
                                           fontsize=ticklabels)

if movie:
    fig.text(right_x,
             bottom_y*0.5,
             f't\n0',
             fontsize=biglabels,
             color='grey',
             ha='right')

# Assign axs objects to variables
# (AxesImage)

artistsgames = np.empty_like(axgames)
artiststraits = np.empty_like(axtraits)
dummy_Zg = np.empty((ext, ext, 4), dtype=np.float32)
dummy_Z = np.empty((nr, nc), dtype=np.float32)
frames = ts
frame0 = ts[-1]

for f, folder in enumerate(folders):
    for a, alpha in enumerate(alphas):
        for l, loges in enumerate(logess):
            ax = axgames[f, a, l]
            artistsgames[f, a, l] = ax.imshow(dummy_Zg,
                                              extent=extent)
    for c in range(len(titles) - 1):
        artiststraits[f, c] = axtraits[f, c].imshow(dummy_Z,
                                                    vmin=0,
                                                    vmax=vmaxs[c])

# Add data and save figure

init(artistsgames,)

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(artiststraits,),
                        blit=True)
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    update(frame0, artiststraits,)
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
