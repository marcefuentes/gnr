#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib import cm
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
titles = ['Games',
          'Severity of\nsocial dilemma', 
          'Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner']
folders = ['given100', 'given95', 'given50']
subfolders = ['p', 'r']

plotsize = 6

# Data

dfprivates = np.empty(len(folders), dtype=object)
dfsocials = np.empty(len(folders), dtype=object)
for f, folder in enumerate(folders):
    filelist = glob(os.path.join(folder, 'none', '*.csv'))
    dfprivates[f] = my.read_files(filelist, False)
    filelist = glob(os.path.join(folder, 'given0', '*.csv'))
    dfsocials[f] = my.read_files(filelist, False)

dftraits = np.empty((len(folders), len(subfolders)), dtype=object)
for f, folder in enumerate(folders):
    for c, subfolder in enumerate(subfolders):
        filelist = glob(os.path.join(folder, subfolder, '*.csv'))
        dftraits[f, c] = my.read_files(filelist, False)

df = dftraits[0, 0]
t = df.Time.max()
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1. - 1./pow(2., logess)
RR, AA = np.meshgrid(rhos, alphas)

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(folders)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabel = plotsize*7
midlabel = plotsize*6
letterlabel = plotsize*5
ticklabel = plotsize*4
xticks = [0, nc/2 - 0.5, nc - 1]
yticks = [0, nr/2 - 0.5, nr - 1]
xmin = logess[0]
xmax = logess[-1]
ymin = alphas[-1]
ymax = alphas[0]
xticklabels = [f'{xmin:.0f}',
               f'{(xmin + xmax)/2.:.0f}',
               f'{xmax:.0f}']
yticklabels = [f'{ymax:.1f}',
               f'{(ymin + ymax)/2:.1f}',
               f'{ymin:.1f}']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig, axs = plt.subplots(nrows=len(folders),
                        ncols=len(titles),
                        figsize=(width, height))


# Apply properties

left_x = axs[0, 0].get_position().x0
right_x = axs[-1, -1].get_position().x1
center_x = (left_x + right_x)/2.
top_y = axs[0, 0].get_position().y1
bottom_y = axs[-1, -1].get_position().y0
center_y = (top_y + bottom_y)/2.
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.3,
              fontsize=biglabel)
fig.supylabel(ylabel,
              x=left_x*0.4,
              y=center_y,
              fontsize=biglabel)

ox = -3/72.
oy = 0/72.
offset = matplotlib.transforms.ScaledTranslation(ox,
                                                 oy,
                                                 fig.dpi_scale_trans)

letter = ord('a')
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
    axs[0, c].set_title(title, pad=plotsize*9, fontsize=midlabel)
    axs[-1, c].set_xticklabels(xticklabels, fontsize=ticklabel)
    for label in axs[-1, c].xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)

# Add data

for f, folder in enumerate(folders):
    lows = my.getZ(t, dftraits[f, 0], 'a2low')
    highs = my.getZ(t, dftraits[f, 0], 'a2high')
    given = dftraits[f, 0].Given.iloc[0]
    T = my.fitness(highs, lows, given, AA, RR)
    R = my.fitness(highs, highs, given, AA, RR)
    P = my.fitness(lows, lows, given, AA, RR)
    S = my.fitness(lows, highs, given, AA, RR)
    Z = my.gamecolors(T, R, P, S)
    axs[f, 0].imshow(Z)

    Zsocial = my.getZ(t, dfsocials[f], 'wmean')
    Z = my.getZ(t, dfprivates[f], 'wmean')
    Z = Zsocial - Z       
    axs[f, 1].imshow(Z)

    for c, trait in enumerate(traits):
        Z = my.getZ(t, dftraits[f, c], trait)
        if 'Grain' in trait:
            Z = 1.0 - Z
        axs[f, c + 2].imshow(Z,
                             vmin=0,
                             vmax=my.a2max) 

# Finish

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
