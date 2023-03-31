#! /usr/bin/env python

from glob import glob
import os
import re
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

traits = ['a2Seen',
          'w']
titles = ['Effort to get $\it{B}$',
          'Fitness']
vmaxs = [my.a2max, my.wmax]
folders = ['none', 'p']

alpha = 0.66
loges = 1.0
movie = False
plotsize = 6
bins = 64

# Add data to figure

def update(t, artists):
    for f, folder in enumerate(folders):
        df = dffrqs[f]
        dfmean = dfmeans[f]
        if movie:
            m = df.Time == t
            df = df.loc[m]
            m = dfmean.Time == t
            dfmean = dfmean.loc[m]
        m = (df.alpha == alpha) & (df.logES == loges)
        d = df.loc[m]
        m = (dfmean.alpha == alpha) & (dfmean.logES == loges)
        dmean = dfmean.loc[m]
        for c, trait in enumerate(traits):
            freq_a = [col for col in d.columns if re.match(fr'^{trait}\d+$', col)]
            y = d.loc[:, freq_a]
            y = y.values[0]
            y = y.flatten()
            artists[f, c].set_ydata(y)
            y = dmean[trait + 'mean'].iloc[0]
            bgcolor = cm.viridis(y/vmaxs[c])
            artists[f, c].axes.set_facecolor(bgcolor)
    if movie:
        fig.texts[1].set_text(f't\n{t}')
    return artists.flatten()

# Data

dffrqs = np.empty(len(folders), dtype=object) 
for f, folder in enumerate(folders):
    filelist = glob(os.path.join(folder, '*.frq'))
    dffrqs[f] = my.read_files(filelist, movie)

dfmeans = np.empty(len(folders), dtype=object) 
for f, folder in enumerate(folders):
    filelist = glob(os.path.join(folder, '*.csv'))
    dfmeans[f] = my.read_files(filelist, movie)

df = dffrqs[1]
ts = df.Time.unique()
rho = 1.0 - 1.0/pow(2, loges)

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(folders)
ylabel = 'Frequency'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*3
xlim = [0, bins - 1]
ylim = [0, 0.2]
xticks = [xlim[0], xlim[1]/2.0, xlim[1]]
yticks = [ylim[0], ylim[1]/2.0, ylim[1]]
xticklabels = np.zeros((len(traits), len(xticks)))
for c, trait in enumerate(traits):
    xticklabels[c] = np.linspace(0.0, vmaxs[c], num=len(xticks))
yticklabels = [f'{ylim[0]:.1f}',
               f'{ylim[1]/2.0:.1f}',
               f'{ylim[1]:.1f}']
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
fig.supylabel(ylabel,
              x=left_x*0.1,
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
    ax.set(xlim=xlim, ylim=ylim)
for f, folder in enumerate(folders):
    axs[f, 0].set_yticklabels(yticklabels, fontsize=ticklabels)
for c, title in enumerate(titles):
    axs[-1, c].set_xlabel(title,
                          labelpad=plotsize*3,
                          fontsize=plotsize*5)
    axs[-1, c].set_xticklabels(xticklabels[c],
                               fontsize=ticklabels)
if movie:
    fig.text(right_x,
             bottom_y*0.5,
             f't\n0',
             fontsize=biglabels,
             color='grey',
             ha='right')

# Assign axs objects to variables
# (Line2D)

artists = np.empty_like(axs) 
x = np.arange(64)
dummy_y = np.zeros_like(x)
frames = ts
frame0 = ts[-1]

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        ax = axs[f, c] 
        artists[f, c], = ax.plot(x,
                                 dummy_y,
                                 c='white',
                                 linewidth=4)

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
