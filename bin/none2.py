#! /usr/bin/env python

from glob import glob
import os
import time

import matplotlib.pyplot as plt
import matplotlib.transforms
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

# Options

titles = ['Effort to get $\it{B}$\nin theory',
          'Effort to get $\it{B}$\nin simulations',
          'Fitness\nin theory',
          'Fitness\nin simulations']
givens_theory = [1.0, 0.95, 0.5, 0.]
givens = ['given100', 'given95', 'given50', 'given00']
vmaxs = [my.a2max,
         my.a2max,
         my.wmax,
         my.wmax]
plotsize = 4

# Data

dfs = np.empty(len(givens), dtype=object)
for g, given in enumerate(givens):
    filelist = glob(os.path.join(given, 'none', '*.csv'))
    dfs[g] = my.read_files(filelist, False)

df = dfs[0]
t = df.Time.max()
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1. - 1./pow(2., logess)
RR, AA = np.meshgrid(rhos, alphas)

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(givens)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabel = plotsize*6
letterlabel = plotsize*5
ticklabel = plotsize*4
xticks = [0, nc/2 - 0.5, nc - 1]
yticks = [0, nr/2 - 0.5, nr - 1]
xmin = my.logesmin
xmax = my.logesmax
ymin = my.alphamin
ymax = my.alphamax
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
              y=bottom_y*0.4,
              fontsize=biglabel)
fig.supylabel(ylabel,
              x=left_x*0.35,
              y=center_y,
              fontsize=biglabel)

ox = 0/72.
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
for g, given in enumerate(givens):
    axs[g, 0].set_yticklabels(yticklabels, fontsize=ticklabel)
for c, title in enumerate(titles):
    axs[0, c].set_title(title, pad=plotsize*10, fontsize=letterlabel)
    axs[-1, c].set_xticklabels(xticklabels,
                               fontsize=ticklabel)
    for label in axs[-1, c].xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)

# Add data to figure

for g, given in enumerate(givens_theory):
    a2_theory = my.a2eq(given, AA, RR)
    axs[g, 0].imshow(a2_theory, vmin=0, vmax=vmaxs[0])
    a2 = my.getZ(t, dfs[g], 'a2Seenmean')
    axs[g, 1].imshow(a2, vmin=0, vmax=vmaxs[1])
    w_theory = my.fitness(a2_theory, a2_theory, given, AA, RR)
    axs[g, 2].imshow(w_theory, vmin=0, vmax=vmaxs[2])
    w = my.getZ(t, dfs[g], 'wmean')
    axs[g, 3].imshow(w, vmin=0, vmax=vmaxs[3])
        
plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')