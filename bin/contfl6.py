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
this_file = os.path.basename(__file__)
file_name = this_file.split('.')[0]

# Options

traits = ['ChooseGrainmean',
          'MimicGrainmean']
titles = ['Sensitivity for\nchoosing partner',
          'Sensitivity for\nmimicking partner']
folders = ['given100', 'given095', 'given050']
subfolders = ['p', 'r']

numaB = 64
theory = False
plotsize = 6
r = 0.1

# Data

filelist = glob(os.path.join('given000', 'none', '*.csv'))
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
xaxis = np.linspace(0.01, my.aBmax - 0.01, num=numaB)

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(folders)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Influence of $\it{B}$'
biglabel = plotsize*7
midlabel = plotsize*6
letterlabel = plotsize*5
ticklabel = plotsize*4
xlim = [0., my.aBmax]
ylim = [0., 1.]
step = int(nr/2)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=len(folders),
                             ncols=len(titles))
axs = np.empty((len(folders),
                len(titles),
                nr,
                nc),
               dtype=object)

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        grid = outergrid[f, c].subgridspec(nrows=nr,
                                           ncols=nc,
                                           wspace=0,
                                           hspace=0)
        axs[f, c] = grid.subplots()

left_x = axs[0, 0, 0, 0].get_position().x0
right_x = axs[-1, -1, -1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0, 0, 0].get_position().y1
bottom_y = axs[-1, -1, -1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.3,
              fontsize=biglabel)
fig.supylabel(ylabel,
              x=left_x*0.1,
              y=center_y,
              fontsize=biglabel)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    ax.set(xlim=xlim, ylim=ylim)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.1)

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        letter = ord('a') + f*len(titles) + c
        axs[f, c, 0, 0].set_title(chr(letter),
                                  fontsize=letterlabel,
                                  pad = 11,
                                  weight='bold',
                                  loc='left')
        if f == 0:
            axs[0, c, 0, 10].set_title(title,
                                       pad=plotsize*9,
                                       fontsize=midlabel)
        for a in range(0, nr, step):
            axs[f, c, a, 0].set(yticks=[ylim[1]/2.0], yticklabels=[])
            if c == 0:
                axs[f, 0, a, 0].set_yticklabels([alphas[a]],
                                                fontsize=ticklabel)
        for e in range(0, nc, step):
            axs[f, c, -1, e].set(xticks=[xlim[1]/2.0], xticklabels=[])
        if folder == folders[-1]:
            for e in range(0, nc, step):
                axs[f, c, -1, e].set_xticklabels([f'{logess[e]:.0f}'],
                                                 fontsize=ticklabel)

# Add data

if theory:
    aBsocials = my.aBeq(0., AA, RR)
else:
    aBsocials = my.getZ(t, dfsocial, 'a2Seenmean')
wsocials = my.getZ(t, dfsocial, 'wmean')

for f, folder in enumerate(folders):

    given = dfprivates[f].Given.iloc[0]
    if theory:
        aBprivates = my.aBeq(given, AA, RR)
    else:
        aBprivates = my.getZ(t, dfprivates[f], 'a2Seenmean')

    wprivates = my.fitness(aBprivates, aBprivates, given, AA, RR)

    Z = np.zeros((len(traits), len(alphas), len(rhos)))
    for c, trait in enumerate(traits):
        Z[c] = my.getZ(t, dftraits[f, c], trait)
        if 'Grain' in trait:
            Z[c] = 1. - Z[c]

    for a, alpha in enumerate(alphas):
        for e, rho in enumerate(rhos):

            aBs = np.full(xaxis.shape, aBprivates[a, e])

            ax = axs[f, 0, a, e]
            ii = my.fitness(xaxis, xaxis, given, alpha, rho)
            ji = my.fitness(aBs, xaxis, given, alpha, rho)
            jj = my.fitness(aBs, aBs, given, alpha, rho)
            ij = my.fitness(xaxis, aBs, given, alpha, rho)
            y = (jj - ii*r - ji*(1. - r))/((1. - r)*(ii - ij - ji + jj))
            mask = xaxis < aBprivates[a, e]
            y[mask] = np.nan
            ax.plot(xaxis, y, color='white', linewidth=0.7)
            #ax.plot(aBprivates[a, e],
            #        wprivates[a, e],
            #        marker='o',
            #        markersize=1.5,
            #        color='white')
            color = cm.viridis(Z[0, a, e]/my.aBmax)
            ax.set_facecolor(color)

            ax = axs[f, 1, a, e]
            y = (my.repeats*my.cost + jj - ji)/(my.repeats*ii - ij - ji + 2.*jj - my.repeats*jj) 
            ax.plot(xaxis, y, color='white', linewidth=0.7)
            ax.plot(aBprivates[a, e],
                    wprivates[a, e],
                    marker='o',
                    markersize=1.5,
                    color='white')
            color = cm.viridis(Z[1, a, e]/my.aBmax)
            ax.set_facecolor(color)

# Finish

plt.savefig(file_name + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
