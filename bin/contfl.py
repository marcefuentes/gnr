#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib.animation import FuncAnimation
from matplotlib import cm
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

predictors = ['Fitness of\npartner choosers',
              'Fitness of\nreciprocators']
traits = ['ChooseGrainmean',
          'MimicGrainmean']
titles_traits = ['Sensitivity for\nchoosing partner',
                 'Sensitivity for\nmimicking partner']
vmaxs = [my.a2max, my.a2max]
folders = ['given100', 'given95', 'given50']
subfolders = ['p', 'r']

numa2 = 64
theory = False
movie = False
plotsize = 6
bgcolor = [0.0, 0.5, 0.7, 1.0]

# Add data to figure

def update(t, lines, images):
    for f, folder in enumerate(folders):
        given = dfprivates[f].Given.iloc[0]
        if theory:
            a2privates = my.a2eq(given, AA, RR)
            ws = my.fitness(a2privates, a2privates, given, AA, RR)
        else: 
            a2privates = my.getZ(t, dfprivates[f], 'a2Seenmean')
            ws = my.getZ(t, dfprivates[f], 'wmean')
        for a, alpha in enumerate(alphas):
            for l, rho in enumerate(rhos):
                w = ws[a, l]
                a2s = np.full(xaxis.shape, a2privates[a, l])
                y = my.fitness(xaxis, a2s, given, alpha, rho)
                color = cm.viridis((my.wmax - y[0])/my.wmax)
                lines[f, 0, a, l].set_ydata(y)
                lines[f, 0, a, l].axes.set_facecolor(color)

                y = my.fitness(xaxis, xaxis, given, alpha, rho) 
                cmap = ListedColormap(['blue', 'yellow'])
                norm = BoundaryNorm([0.0, w], cmap.N)
                points = np.array([xaxis, y]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                lc = LineCollection(segments,
                                    cmap=cmap,
                                    norm=norm,
                                    linewidth=2,
                                    #facecolor= 'white',
                                    array=y)
                ax = lines[f, 1, a, l].axes
                ax.add_collection(lc)
        for c, trait in enumerate(traits):
            Z = my.getZ(t, dftraits[f, c], trait)
            if 'Grain' in trait:
                Z = 1.0 - Z
            images[f, c].set_array(Z)
    if movie:
        fig.texts[2].set_text(f't\n{t}')
    return np.concatenate([lines.flatten(), images.flatten()]) 

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
xaxis = np.linspace(0.01, my.a2max - 0.01, num=numa2)

# Figure properties

width = plotsize*(len(predictors) + len(traits))
height = plotsize*len(folders)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
xlim = [0.0, my.a2max]
ylim = [0.0, my.wmax]
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

left_x = axlines[0, 0, 0, 0].get_position().x0
right_x = aximages[-1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = aximages[0, 0].get_position().y1
bottom_y = aximages[-1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.3,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=left_x*0.5,
              y=center_y,
              fontsize=biglabels)

ox = -3/72.; oy = 0/72.
offset = matplotlib.transforms.ScaledTranslation(ox,
                                                 oy,
                                                 fig.dpi_scale_trans)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])

letter = ord('a')
letterposition = 1.035
for f, folder in enumerate(folders):
    for c, predictor in enumerate(predictors): 
        for a, alpha in enumerate(alphas):
            for l, loges in enumerate(logess):
                for axis in ['top','bottom','left','right']:
                    axlines[f, c, a, l].spines[axis].set_linewidth(0.1)
                axlines[f, c, a, l].set(xlim=xlim, ylim=ylim)
                #axlines[f, c, a, l].set_facecolor(bgcolor)
        axlines[f, c, 0, 0].set_title(chr(letter),
                                      pad=plotsize*5/3,
                                      fontsize=plotsize*5,
                                      weight='bold')
        letter += 1
        if f == 0:
            axlines[0, c, 0, 10].set_title(predictor,
                                           pad=plotsize*9,
                                           fontsize=plotsize*5)
        if c == 0:
            for a in range(0, nr, step):
                axlines[f, 0, a, 0].set(yticks=[my.wmax/2.0])
                axlines[f, 0, a, 0].set_yticklabels([f'{alphas[a]:.1f}'],
                                                    rotation='horizontal',
                                                    horizontalalignment='right',
                                                    verticalalignment='center',
                                                    fontsize=ticklabels)
        for l in range(0, nc, step):
            axlines[f, c, -1, l].set(xticks=[my.a2max/2.0], xticklabels=[]) 
        if folder == folders[-1]:
            for l in range(0, nc, step):
                axlines[-1, c, -1, l].set_xticklabels([f'{logess[l]:.0f}'],
                                                      fontsize=ticklabels)
    for c, titles_trait in enumerate(titles_traits):
        aximages[f, c].text(0,
                            letterposition,
                            chr(letter),
                            transform=aximages[f, c].transAxes,
                            fontsize=plotsize*5,
                            weight='bold')
        letter += 1
        aximages[f, c].set(xticks=xticks, yticks=yticks)
        aximages[f, c].set(xticklabels=[], yticklabels=[])
        if f == 0:
            aximages[0, c].set_title(titles_trait,
                                     pad=plotsize*9,
                                     fontsize=plotsize*5)
        if folder == folders[-1]:
            aximages[-1, c].set_xticklabels(xticklabels,
                                            fontsize=ticklabels)
            for label in aximages[-1, c].xaxis.get_majorticklabels():
                label.set_transform(label.get_transform() + offset)
if movie:
    fig.text(right_x,
             bottom_y*0.5,
             f't\n0',
             fontsize=biglabels,
             color='grey',
             ha='right')

# Assign axs objects to variables
# (Line2D)
# (AxesImage)

lines = np.empty_like(axlines)
images = np.empty_like(aximages)
dummy_Z = np.empty((nr, nc), dtype=float)
dummy_y = np.empty_like(xaxis)
frames = ts
frame0 = ts[-1]

for f, folder in enumerate(folders):
    for c, predictor in enumerate(predictors):
        for a, alpha in enumerate(alphas):
            for l, loges in enumerate(logess):
                ax = axlines[f, c, a, l]
                lines[f, c, a, l], = ax.plot(xaxis,
                                             dummy_y,
                                             color='white')
    for c, trait in enumerate(traits):
        ax = aximages[f, c]
        images[f, c] = ax.imshow(dummy_Z,
                                 vmin=0,
                                 vmax=vmaxs[c])

# Add data and save figure

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(lines, images,),
                        blit=True)
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    update(frame0, lines, images,)
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
