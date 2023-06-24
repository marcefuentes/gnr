#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib import cm
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
this_file = os.path.basename(__file__)
file_name = this_file.split(".")[0]

# Options

traits = ["ChooseGrainmean",
          "MimicGrainmean"]
titles = ["Games",
          "Sensitivity for\nchoosing partner",
          "Sensitivity for\nmimicking partner"]
folders = ["given100", "given095", "given050"]
subfolders = ["p", "r"]

movie = False
plotsize = 6

# Add data to figure

def init(artists):

    for f, folder in enumerate(folders):
        lows = my.getZ(ts[0], dftraits[f, 0], "aBlow")
        highs = my.getZ(ts[0], dftraits[f, 0], "aBhigh")
        given = dftraits[f, 0].Given.iloc[0]
        T = my.fitness(highs, lows, given, AA, RR)
        R = my.fitness(highs, highs, given, AA, RR)
        P = my.fitness(lows, lows, given, AA, RR)
        S = my.fitness(lows, highs, given, AA, RR)
        #Ma = np.maximum.reduce([T, R, P, S])
        #Mi = np.minimum.reduce([T, R, P, S])
        #Tn = (T - Mi)/(Ma - Mi)
        #Rn = (R - Mi)/(Ma - Mi)
        #Pn = (P - Mi)/(Ma - Mi)
        #Sn = (S - Mi)/(Ma - Mi)
        Tn, Rn, Pn, Sn = T, R, P, S
        y = np.stack((Tn, Rn, Pn, Sn), axis=-1)
        linecolor = np.full(highs.shape, "white")
        red = np.full(highs.shape, "red")
        m = lows > highs
        linecolor[m] = red[m]

        Zg = my.gamecolors(T, R, P, S)
        for c, title in enumerate(titles):
            for (a, l, i), _ in np.ndenumerate(y):
                artists[f, c, a, l].set_ydata(y[a, l])
                if c == 0:
                    artists[f, c, a, l].axes.set_facecolor(Zg[a, l])
                lcolor = linecolor[a, l] 
                artists[f, c, a, l].set_color(lcolor)
                artists[f, c, a, l].set_markerfacecolor(lcolor)

    return artists.flatten()

def update(t, artists):
    for f, folder in enumerate(folders):
        for c, trait in enumerate(traits):
            Z = my.getZ(t, dftraits[f, c], trait)
            if "Grain" in trait:
                Z = 1.0 - Z
            for (a, l), _ in np.ndenumerate(Z):
                bgcolor = cm.viridis(Z[a, l]/my.aBmax)
                artists[f, c + 1, a, l].axes.set_facecolor(bgcolor)
    if movie:
        fig.texts[2].set_text(f"t\n{t}")
    return artists.flatten()

# Data

dftraits = np.empty((len(folders), len(subfolders)), dtype=object)
for f, folder in enumerate(folders):
    for c, subfolder in enumerate(subfolders):
        filelist = glob(os.path.join(folder, subfolder, "*.csv"))
        dftraits[f, c] = my.read_files(filelist, movie)

df = dftraits[0, 0]
ts = df.Time.unique()
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1.0 - 1.0/pow(2.0, logess)
RR, AA = np.meshgrid(rhos, alphas)

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(folders)
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
biglabel = plotsize*7
midlabel = plotsize*6
letterlabel = plotsize*5
ticklabel = plotsize*4
xlim = [0, 5]
#ylim = [-0.1, 1.1]
ylim = [0.0, 2.0]
step = int(nr/2)
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

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
              x=left_x*0.4,
              y=center_y,
              fontsize=biglabel)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    ax.set(xlim=xlim, ylim=ylim)
    for axis in ["top","bottom","left","right"]:
        ax.spines[axis].set_linewidth(0.1)

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        letter = ord("a") + f*len(titles) + c
        axs[f, c, 0, 0].set_title(chr(letter),
                                  fontsize=lettersize,
                                  pad = 11,
                                  weight="bold",
                                  loc="left")
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
                axs[f, c, -1, e].set_xticklabels([f"{logess[e]:.0f}"],
                                                 fontsize=ticklabel)
if movie:
    fig.text(right_x,
             bottom_y*0.5,
             "t\n0",
             fontsize=biglabels,
             color="grey",
             ha="right")

# Assign axs objects to variables
# (Line2D)

artists = np.empty_like(axs) 
xaxis = [1, 2, 3, 4]
dummy_y = np.zeros_like(xaxis)
frames = ts
frame0 = ts[-1]

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        for a, alpha in enumerate(alphas):
            for l, loges in enumerate(logess):
                ax = axs[f, c, a, l] 
                artists[f, c, a, l], = ax.plot(xaxis,
                                               dummy_y,
                                               linewidth=1,
                                               marker="o",
                                               markersize=plotsize/3)

# Add data and save figure

init(artists,)

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(artists,),
                        blit=True)
    ani.save(file_name + ".mp4", writer="ffmpeg", fps=10)
else:
    update(frame0, artists,)
    plt.savefig(file_name + ".png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
