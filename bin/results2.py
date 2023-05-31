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
filename = this_file.split('.')[0]

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
folders = ['given100', 'given095', 'given050', 'given000']
subfolder = 'r'

movie = False
plotsize = 4

# Add data to figure

def update(t, artists):
    wsocial = my.getZ(t, dfs[-1], 'wmean')
    for f, folder in enumerate(folders):
        for c, trait in enumerate(traits):
            Z = my.getZ(t, dfs[f], trait)
            if 'Grain' in trait:
                Z = 1. - Z
            if 'gain' in titles[c]:
                wnull = my.getZ(t, dfnulls[f], 'wmean')
                Z = Z - wnull
            if 'deficit' in titles[c]:
                Z = wsocial - Z
            artists[f, c].set_array(Z) 
    if movie:
        fig.texts[2].set_text(f't\n{t}')
    return artists.flatten()

# Data

dfnulls = np.empty(len(folders), dtype=object) 
for f, folder in enumerate(folders):
    filelist = glob(os.path.join(folder, 'none', '*.csv'))
    dfnulls[f] = my.read_files(filelist, movie)

dfs = np.empty(len(folders), dtype=object) 
for f, folder in enumerate(folders):
    filelist = glob(os.path.join(folder, subfolder, '*.csv'))
    dfs[f] = my.read_files(filelist, movie)

df = dfs[1]
ts = df.Time.unique()
nr = df['alpha'].nunique()
nc = df['logES'].nunique()

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(folders)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Influence of $\it{B}$'
biglabel = plotsize*7
letterlabel = plotsize*5
ticklabel = plotsize*4
xticks = [0, nc/2 - 0.5, nc - 1]
yticks = [0, nr/2 - 0.5, nr - 1]
xmin = df['logES'].min()
xmax = df['logES'].max()
ymin = df['alpha'].min()
ymax = df['alpha'].max()
xticklabels = [f'{xmin:.0f}',
               f'{(xmin + xmax)/2.:.0f}',
               f'{xmax:.0f}']
yticklabels = [f'{ymax:.1f}',
               f'{(ymin + ymax)/2.:.1f}',
               f'{ymin:.1f}']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig, axs = plt.subplots(nrows=len(folders),
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
              y=bottom_y*0.3,
              fontsize=biglabel)
fig.supylabel(ylabel,
              x=left_x*0.3,
              y=center_y,
              fontsize=biglabel)

ox = -2/72.
oy = 0/72.
offset = matplotlib.transforms.ScaledTranslation(ox, oy, fig.dpi_scale_trans)

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
for f, folder in enumerate(folders):
    axs[f, 0].set_yticklabels(yticklabels, fontsize=ticklabel)
for c, title in enumerate(titles):
    axs[0, c].set_title(title, pad=plotsize*10, fontsize=letterlabel)
    axs[-1, c].set_xticklabels(xticklabels,
                               fontsize=ticklabel)
    for label in axs[-1, c].xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)

if movie:
    fig.text(right_x,
             bottom_y*0.5,
             't\n0',
             fontsize=biglabel,
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
