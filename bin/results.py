#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import matplotlib.transforms
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
this_file = os.path.basename(__file__)
file_name = this_file.split('.')[0]

# Options

traits = ['ChooseGrainmean',
          'MimicGrainmean',
          'wmean',
          'wmean']
titles = ['Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner',
          'Fitness gain',
          'Fitness deficit']
vmaxs = [my.aBmax,
         my.aBmax,
         my.wmax,
         my.wmax]
givens = ['given100', 'given095', 'given050', 'given000']
folders = ['none', 'p', 'p8', 'p8r', 'pr', 'r']

movie = False
plotsize = 4

# Add data to figure

def update(t, artists):
    wsocial = my.getZ(t, dfs[-1], 'wmean')
    for g, given in enumerate(givens):
        for c, trait in enumerate(traits):
            Z = my.getZ(t, dfs[g], trait)
            if 'Grain' in trait:
                Z = 1. - Z
            if 'gain' in titles[c]:
                wnull = my.getZ(t, dfnulls[g], 'wmean')
                Z = Z - wnull
            if 'deficit' in titles[c]:
                Z = wsocial - Z
            artists[g, c].set_array(Z) 
    if movie:
        fig.texts[2].set_text(f't\n{t}')
    return artists.flatten()

# Data without partner choice or reciprocity

dfnulls = np.empty(len(givens), dtype=object) 
for g, given in enumerate(givens):
    filelist = glob(os.path.join('none', given, '*.csv'))
    dfnulls[g] = my.read_files(filelist, movie)

df = dfnulls[0]
ts = df.Time.unique()
nr = df.alpha.nunique()
nc = df.logES.nunique()

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(givens)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Influence of $\it{B}$'
biglabel = plotsize*7
letterlabel = plotsize*6
ticklabel = plotsize*5
xticks = [0, nc/2 - 0.5, nc - 1]
yticks = [0, nr/2 - 0.5, nr - 1]
xmin = df.logES.min()
xmax = df.logES.max()
ymin = df.alpha.min()
ymax = df.alpha.max()
xticklabels = [f'{xmin:.0f}',
               f'{(xmin + xmax)/2.:.0f}',
               f'{xmax:.0f}']
yticklabels = [f'{ymax:.1f}',
               f'{(ymin + ymax)/2.:.1f}',
               f'{ymin:.1f}']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig, axs = plt.subplots(nrows=len(givens),
                        ncols=len(titles),
                        figsize=(width, height))

left_x = axs[0, 0].get_position().x0
right_x = axs[-1, -1].get_position().x1
center_x = (left_x + right_x) / 2.
top_y = axs[0, 0].get_position().y1
bottom_y = axs[-1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2.
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y - 1.2/height,
              fontsize=biglabel)
fig.supylabel(ylabel,
              x=left_x - 1.45/width,
              y=center_y,
              fontsize=biglabel)

letterposition = 1.035
for i, ax in enumerate(fig.get_axes()):
    ax.set(xticks=xticks, yticks=yticks)
    ax.set(xticklabels=[], yticklabels=[])
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.1)
    letter = ord('a') + i
    ax.text(0,
            letterposition,
            chr(letter),
            transform=ax.transAxes,
            fontsize=letterlabel,
            weight='bold')
for g, given in enumerate(givens):
    axs[g, 0].set_yticklabels(yticklabels, fontsize=ticklabel)
for c, title in enumerate(titles):
    axs[0, c].set_title(title, pad=plotsize*10, fontsize=letterlabel)
    axs[-1, c].set_xticklabels(xticklabels,
                               fontsize=ticklabel)
if movie:
    fig.text(right_x,
             bottom_y*0.5,
             't\n0',
             fontsize=biglabel,
             color='grey',
             ha='right')

# Assign axs objects to variables
# (AxesImage)

dfs = np.empty(len(givens), dtype=object) 
artists = np.empty_like(axs) 
dummy_Z = np.zeros((nr, nc))
frames = ts
frame0 = ts[-1]

for folder in folders:
    for g, given in enumerate(givens):
        filelist = glob(os.path.join(folder, given, '*.csv'))
        dfs[g] = my.read_files(filelist, movie)
        for c, title in enumerate(titles):
            artists[g, c] = axs[g, c].imshow(dummy_Z,
                                             vmin=0,
                                             vmax=vmaxs[c])

    # Add data and save figure

    if movie:
        ani = FuncAnimation(fig,
                            update,
                            frames=frames,
                            fargs=(artists,),
                            blit=True)
        ani.save(file_name + '_' + folder + '.mp4', writer='ffmpeg', fps=10)
    else:
        update(frame0, artists,)
        plt.savefig(file_name + '_' + folder + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
