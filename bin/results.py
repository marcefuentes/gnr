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

traits = ['a2Seenmean',
          'ChooseGrainmean',
          'MimicGrainmean',
          'wmean']
titles = ['Effort to get $\it{B}$',
          'Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner',
          'Fitness']
vmaxs = [my.a2max, my.a2max, my.a2max, 2.0]
folders = ['given0', 'none', 'p', 'r', 'pr', 'p8r']

movie = False
plotsize = 4

# Add data to figure

def update(t, images):
    for f, folder in enumerate(folders):
        df = dfs[f]
        if movie:
            m = df.Time == t
            df = df.loc[m]
        for c, trait in enumerate(traits):
            Z = pd.pivot_table(df,
                               values=trait,
                               index='alpha',
                               columns='logES')
            Z = Z.sort_index(axis=0, ascending=False)
            if 'Grain' in trait:
                Z = 1.0 - Z
            images[f, c].set_array(Z) 
    if movie:
        fig.texts[2].set_text(f't\n{t}')
    return images.flatten()

# Get data

dfs = np.empty(len(folders), dtype=object) 
for i, folder in enumerate(folders):
    filelist = glob(os.path.join(folder, '*.csv'))
    df = [my.read_file(file, movie) for file in filelist]
    dfs[i] = pd.concat(df, ignore_index=True)

df = dfs[1]
ts = df.Time.unique()
nr = df['alpha'].nunique()
nc = df['logES'].nunique()

# Figure properties

width = plotsize*len(traits)
height = plotsize*len(folders)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
xticks = [0, nc/2-0.5, nc-1]
yticks = [0, nr/2-0.5, nr-1]
xmin = df['logES'].min()
xmax = df['logES'].max()
ymin = df['alpha'].min()
ymax = df['alpha'].max()
xticklabels = [f'{xmin:2.0f}',
               f'{(xmin + xmax)/2.0:2.0f}',
               f'{xmax:2.0f}']
yticklabels = [f'{ymax:3.1f}',
               f'{(ymin + ymax)/2.0:3.1f}',
               f'{ymin:3.1f}']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig, axs = plt.subplots(nrows=len(folders),
                        ncols=len(traits),
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
    axs[-1, c].set_xticklabels(xticklabels, fontsize=ticklabels)

if movie:
    fig.text(right_x,
             bottom_y*0.5,
             f't\n0',
             fontsize=biglabels,
             color='grey',
             ha='right')

# Assign axs objects to variables
# (AxesImage artists to images)

images = np.empty(axs.shape, dtype=object) 
dummy_Z = np.empty((nr, nc), dtype=np.float32)

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        images[f, c] = axs[f, c].imshow(dummy_Z,
                                        vmin=0,
                                        vmax=vmaxs[c])

# Add data and save figure

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=ts,
                        fargs=(images,),
                        blit=True)
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    update(ts[-1], images,)
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
