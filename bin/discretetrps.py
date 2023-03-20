#! /usr/bin/env python

from glob import glob
import os
import time

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
traitvmaxs = [my.a2max,
              my.a2max]
folders = ['0', '1']
subfolders = ['p', 'r']

plotsize = 4

# Add data to figure

def figdata(images, lines):
    for i, folder in enumerate(folders):
        df = dfs[i, 0]
        given = df.Given.iloc[0]
        lows = pd.pivot_table(df,
                              values='a2low',
                              index=['alpha'],
                              columns=['logES'])
        lows = lows.sort_index(axis=0, ascending=False)
        lows = lows.to_numpy()
        highs = pd.pivot_table(df,
                               values='a2high',
                               index=['alpha'],
                               columns=['logES'])
        highs = highs.sort_index(axis=0, ascending=False)
        highs = highs.to_numpy()
        T = my.fitness(highs, lows, given, AA, RR)
        R = my.fitness(highs, highs, given, AA, RR)
        P = my.fitness(lows, lows, given, AA, RR)
        S = my.fitness(lows, highs, given, AA, RR)

        Z = my.gamecolors(T, R, P, S)
        images[i, 0].set_array(Z)
        for j, trait in enumerate(traits):
            df = dfs[i, j]
            Z = pd.pivot_table(df,
                               values=trait,
                               index=['alpha'],
                               columns=['logES'])
            Z = Z.sort_index(axis=0, ascending=False)
            if 'Grain' in trait:
                Z = 1.0 - Z
            images[i, j+1].set_array(Z) 
    return images.flatten()

# Get data

def read_file(file):
    df = pd.read_csv(file)
    return df.tail(1)

dfs = np.empty((len(folders), len(subfolders)), dtype=object) 
for i, folder in enumerate(folders):
    for j, subfolder in enumerate(subfolders):
        filelist = glob(os.path.join(folder, subfolder, '*.csv'))
        d = list(map(read_file, filelist))
        dfs[i, j] = pd.concat(d, ignore_index=True)

df = dfs[0, 0]
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
xticks = [0, nc/2-0.5, nc-1]
yticks = [0, nr/2-0.5, nr-1]
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
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig, axs = plt.subplots(nrows=len(folders),
                        ncols=len(titles),
                        figsize=(width, height))
images = np.empty(axs.shape, dtype=object) 
dummy_Zg = np.zeros((nr, nc, 4), dtype=np.float32)
dummy_Z = np.empty((nr, nc), dtype=np.float32)

left_x = axs[0, 0].get_position().x0
right_x = axs[-1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0].get_position().y1
bottom_y = axs[-1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.5,
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
for i, folder in enumerate(folders):
    axs[i, 0].set_yticklabels(yticklabels, fontsize=ticklabels)
for j, title in enumerate(titles):
    axs[0, j].set_title(title, pad=plotsize*9, fontsize=plotsize*5)
    axs[-1, j].set_xticklabels(xticklabels, fontsize=ticklabels)

# Assign axs objects to variables
# (Line2D objects to lines)
# (AxesImage objects to "images")

for i, folder in enumerate(folders):

    images[i] = axs[i, 0].imshow(dummy_Zg)

    for j, trait in enumerate(traits):
        images[i, j + 1] = axs[i, j + 1].imshow(dummy_Z,
                                                vmin=0,
                                                vmax=traitvmaxs[j])

# Add data and save figure

figdata(images,)
plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')