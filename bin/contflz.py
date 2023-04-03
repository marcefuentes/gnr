#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.pyplot as plt
import matplotlib.transforms
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

# Options

traits = ['ChooseGrainmean',
          'MimicGrainmean']
titles = ['Fitness',
          'Fitness',
          'Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner']
vmaxs = [my.wmax/2.0, my.wmax/2.0, my.a2max, my.a2max]
folders = ['given100', 'given95', 'given50']
subfolders = ['p', 'r']

dx = 0.01
theory = False
movie = False
plotsize = 6

# Add data to figure

def update(t, artists):
    for f, folder in enumerate(folders):
        given = dfprivates[f].Given.iloc[0]
        if theory:
            a2privates = my.a2eq(given, AA, RR)
            ws = my.fitness(a2privates, a2privates, given, AA, RR)
        else: 
            a2privates = my.getZ(t, dfprivates[f], 'a2Seenmean')
            #ws = my.getZ(t, dfprivates[f], 'wmean')
            ws = my.fitness(a2privates, a2privates, given, AA, RR)

        Z = my.wmax/2.0 - my.fitness(dx, a2privates, given, AA, RR)
        m = Z < 0.0
        Z[m] = zeros[m]
        artists[f, 0].set_array(Z)

        Z = my.fitness(a2privates, a2privates, given, AA, RR)
        Z = Z - my.fitness(dx, dx, given, AA, RR)
        m = Z < 0.0
        Z[m] = zeros[m]
        artists[f, 1].set_array(Z)
        for c, trait in enumerate(traits):
            Z = my.getZ(t, dftraits[f, c], trait)
            if 'Grain' in trait:
                Z = 1.0 - Z
            artists[f, c + 2].set_array(Z)
    if movie:
        fig.texts[2].set_text(f't\n{t}')
    return artists.flatten() 

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
zeros = np.zeros(AA.shape)
dx = np.full(AA.shape, dx)

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(folders)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
xticks = [0, nc/2 - 0.5, nc - 1]
yticks = [0, nr/2 - 0.5, nr - 1]
xmin = logess[0]
xmax = logess[-1]
ymin = alphas[-1]
ymax = alphas[0]
xticklabels = [f'{xmin:.0f}',
               f'{(xmin + xmax)/2.0:.0f}',
               f'{xmax:.0f}']
yticklabels = [f'{ymax:.1f}',
               f'{(ymin + ymax)/2.0:.1f}',
               f'{ymin:.1f}']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig, axs = plt.subplots(nrows=len(folders),
                        ncols=len(titles),
                        figsize=(width, height))

left_x = axs[0, 0].get_position().x0
right_x = axs[-1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0].get_position().y1
bottom_y = axs[-1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.6,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=left_x*0.4,
              y=center_y,
              fontsize=biglabels)

ox = -2/72.; oy = 0/72.
offset = matplotlib.transforms.ScaledTranslation(ox, oy, fig.dpi_scale_trans)

letterposition = 1.035
for i, ax in enumerate(fig.get_axes()):
    ax.set(xticks=xticks, yticks=yticks)
    ax.set(xticklabels=[], yticklabels=[])
    letter = ord('a') + i
    ax.text(0,
            letterposition,
            chr(letter),
            transform=ax.transAxes,
            fontsize=plotsize*5,
            weight='bold')
for f, folder in enumerate(folders):
    axs[f, 0].set_yticklabels(yticklabels, fontsize=ticklabels)
for c, title in enumerate(titles):
    axs[0, c].set_title(title, pad=plotsize*10, fontsize=plotsize*5)
    axs[-1, c].set_xticklabels(xticklabels,
                               fontsize=ticklabels)
    for label in axs[-1, c].xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)

if movie:
    fig.text(right_x,
             bottom_y*0.5,
             f't\n0',
             fontsize=biglabels,
             color='grey',
             ha='right')

# Assign axs objects to variables
# (AxesImage)

artists = np.empty_like(axs) 
dummy_Z = np.empty((nr, nc), dtype=float)
frames = ts
frame0 = ts[-1]

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        artists[f, c] = axs[f, c].imshow(dummy_Z,
                                         vmin=0,
                                         vmax=vmaxs[c])

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
