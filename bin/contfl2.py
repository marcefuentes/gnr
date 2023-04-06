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

predictors = ['$\it{T}$, $\it{S}$',
              '$\it{R}$, $\it{P}$']
traits = ['ChooseGrainmean',
          'MimicGrainmean']
titles_traits = ['Sensitivity for\nchoosing partner',
                 'Sensitivity for\nmimicking partner']
vmaxs = [my.a2max, my.a2max]
folders = ['given100', 'given95', 'given50']
subfolders = ['p', 'r']

numa2 = 64
theory = False
plotsize = 6

# Data

filelist = glob(os.path.join('given00', 'none', '*.csv'))
dfsocial = my.read_files(filelist, False)

dfprivates = np.empty(len(folders), dtype=object)
for f, folder in enumerate(folders):
    filelist = glob(os.path.join(folder, 'none', '*.csv'))
    dfprivates[f] = my.read_files(filelist, False)

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
xaxis = np.linspace(0.01, my.a2max - 0.01, num=numa2)

# Figure properties

width = plotsize*(len(predictors) + len(traits))
height = plotsize*len(folders)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabel = plotsize*7
midlabel = plotsize*6
letterlabel = plotsize*5
ticklabel = plotsize*4
xlim = [0., my.a2max]
ylim = [0., my.wmax]
step = int(nr/2)
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

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=len(folders),
                             ncols=len(predictors) + len(traits))
axlines = np.empty((len(folders),
                    len(predictors),
                    nr,
                    nc),
                   dtype=object)
aximages = np.empty((len(folders),
                     len(traits)), 
                    dtype=object)

for f, folder in enumerate(folders):
    for c, predictor in enumerate(predictors):
        grid = outergrid[f, c].subgridspec(nrows=nr,
                                           ncols=nc,
                                           wspace=0,
                                           hspace=0)
        axlines[f, c] = grid.subplots()
    for c, trait in enumerate(traits):
        grid = outergrid[f, len(predictors) + c].subgridspec(nrows=1,
                                                             ncols=1)
        aximages[f, c] = grid.subplots()

# Apply properties

left_x = axlines[0, 0, 0, 0].get_position().x0
right_x = aximages[-1, -1].get_position().x1
center_x = (left_x + right_x)/2.
top_y = aximages[0, 0].get_position().y1
bottom_y = aximages[-1, -1].get_position().y0
center_y = (top_y + bottom_y)/2.
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.3,
              fontsize=biglabel)
fig.supylabel(ylabel,
              x=left_x*0.5,
              y=center_y,
              fontsize=biglabel)

ox = -3/72.
oy = 0/72.
offset = matplotlib.transforms.ScaledTranslation(ox,
                                                 oy,
                                                 fig.dpi_scale_trans)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])

letter = ord('a')
letterposition = 1.035
initial_linewidth = axlines[0, 0, 0, 0].spines['top'].get_linewidth()
for f, folder in enumerate(folders):
    for c, predictor in enumerate(predictors):
        for axis in ['top', 'bottom', 'left', 'right']:
            aximages[f, c].spines[axis].set_linewidth(0.1)
        for a, alpha in enumerate(alphas):
            for e, loges in enumerate(logess):
                for axis in ['top', 'bottom', 'left', 'right']:
                    axlines[f, c, a, e].spines[axis].set_linewidth(0.1)
                axlines[f, c, a, e].set(xlim=xlim, ylim=ylim)
        axlines[f, c, 0, 0].set_title(chr(letter),
                                      pad=plotsize*5/3,
                                      fontsize=letterlabel,
                                      weight='bold')
        letter += 1
        if f == 0:
            axlines[0, c, 0, 10].set_title(predictor,
                                           pad=plotsize*9,
                                           fontsize=midlabel)
        if c == 0:
            for a in range(0, nr, step):
                axlines[f, 0, a, 0].set(yticks=[my.wmax/2.])
                axlines[f, 0, a, 0].set_yticklabels([f'{alphas[a]:.1f}'],
                                                    rotation='horizontal',
                                                    horizontalalignment='right',
                                                    verticalalignment='center',
                                                    fontsize=ticklabel)
        for e in range(0, nc, step):
            axlines[f, c, -1, e].set(xticks=[my.a2max/2.], xticklabels=[]) 
        if folder == folders[-1]:
            for e in range(0, nc, step):
                axlines[-1, c, -1, e].set_xticklabels([f'{logess[e]:.0f}'],
                                                      fontsize=ticklabel)
    for c, titles_trait in enumerate(titles_traits):
        aximages[f, c].text(0,
                            letterposition,
                            chr(letter),
                            transform=aximages[f, c].transAxes,
                            fontsize=letterlabel,
                            weight='bold')
        letter += 1
        aximages[f, c].set(xticks=xticks, yticks=yticks)
        aximages[f, c].set(xticklabels=[], yticklabels=[])
        if f == 0:
            aximages[0, c].set_title(titles_trait,
                                     pad=plotsize*9,
                                     fontsize=midlabel)
        if folder == folders[-1]:
            aximages[-1, c].set_xticklabels(xticklabels,
                                            fontsize=ticklabel)
            for label in aximages[-1, c].xaxis.get_majorticklabels():
                label.set_transform(label.get_transform() + offset)

# Add data

for f, folder in enumerate(folders):

    given = dfprivates[f].Given.iloc[0]
    if theory:
        a2privates = my.a2eq(given, AA, RR)
    else:
        a2privates = my.getZ(t, dfprivates[f], 'a2Seenmean')
        wmean = my.getZ(t, dfprivates[f], 'wmean')
    wpp = my.fitness(a2privates, a2privates, given, AA, RR)

    for a, alpha in enumerate(alphas):
        for e, rho in enumerate(rhos):

            w = wpp[a, e]
            a2s = np.full(xaxis.shape, a2privates[a, e])

            ax = axlines[f, 0, a, e]
            y = my.fitness(xaxis, a2s, given, alpha, rho)
            wsum = 0.
            j = 0
            for i in range(numa2):
                if xaxis[i] < a2s[i]:
                    wsum += w - y[i]
                    j += 1
            if wsum < 0:
                wsum = 0
            color = cm.viridis(wsum/(w*j))
            ax.plot(xaxis, y, color='black')
            ax.set_facecolor(color)

            ax = axlines[f, 1, a, e]
            y = my.fitness(xaxis, xaxis, given, alpha, rho)
            wsum = 0.
            for i in range(numa2):
                wsum += y[i] - w
            if wsum < 0:
                wsum = 0
            color = cm.viridis(wsum/(w*numa2))
            ax.plot(xaxis, y, color='black')
            ax.set_facecolor(color)

    for c, trait in enumerate(traits):
        Z = my.getZ(t, dftraits[f, c], trait)
        if 'Grain' in trait:
            Z = 1. - Z
        aximages[f, c].imshow(Z,
                              vmin=0,
                              vmax=vmaxs[c])

# Finish

plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
