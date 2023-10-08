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
this_file = os.path.basename(__file__)
file_name = this_file.split(".")[0]

# Options

traits = ["a2Seen",
          "ChooseGrain",
          "MimicGrain",
          "w"]
titles = ["Production of $\it{B}$",
          "Sensitivity for\nchoosing partner",
          "Sensitivity for\nmimicking partner",
          "Fitness"]
vmaxs = [my.aBmax,
         my.aBmax,
         my.aBmax,
         my.wmax]
givens = ["given100", "given095", "given050", "given000"]
folders = ["none", "p", "p8", "pr8", "pr", "r"]
folders = ["r"]

movie = False
plotsize = 8
bins = 64

# Add data to figure

def update(t, artists):
    for g, given in enumerate(givens):
        df = dffrqs[g]
        if movie:
            m = df.Time == t
            df = df.loc[m]
        for c, trait in enumerate(traits):
            Z = my.getZ(t, dfmeans[g], trait + "mean")
            if "Grain" in trait:
                Z = 1.0 - Z
            for a, alpha in enumerate(alphas):
                for e, loges in enumerate(logess):
                    m = (df.alpha == alpha) & (df.logES == loges)
                    d = df.loc[m]
                    freq_a = [col for col in d.columns if re.match(fr"^{trait}\d+$", col)]
                    y = d.loc[:, freq_a]
                    y = y.values[0]
                    y = y.flatten()
                    artists[g, c, a, e].set_ydata(y)
                    bgcolor = cm.viridis(Z[a, e]/vmaxs[c])
                    artists[g, c, a, e].axes.set_facecolor(bgcolor)
    if movie:
        fig.texts[2].set_text(f"t\n{t}")
    return artists.flatten()

# Data without partner choice or reciprocity

filelist = glob(os.path.join(folders[0], givens[0], "*.csv"))
df = my.read_files(filelist, movie)
ts = df.Time.unique()
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(givens)
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
biglabel = plotsize*6
letterlabel = plotsize*5
ticklabel = plotsize*4
xlim = [-2, bins + 1]
ylim = [0, 0.25]
step = int(nr/2)
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Create figure

axs = np.empty((len(givens),
                len(titles),
                nr,
                nc),
                dtype=object)

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=len(givens),
                             ncols=len(titles))
for g, given in enumerate(givens):
    for c, title in enumerate(titles):
        grid = outergrid[g, c].subgridspec(nrows=nr,
                                           ncols=nc,
                                           wspace=0,
                                           hspace=0)
        axs[g, c] = grid.subplots()

left_x = axs[0, 0, 0, 0].get_position().x0
right_x = axs[-1, -1, -1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0, 0, 0].get_position().y1
bottom_y = axs[-1, -1, -1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y - 1.4/height,
              fontsize=biglabel)
fig.supylabel(ylabel,
              x=left_x - 1.5/width,
              y=center_y,
              fontsize=biglabel)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    ax.set(xlim=xlim, ylim=ylim)
    for axis in ["top","bottom","left","right"]:
        ax.spines[axis].set_linewidth(0.1)

for g, given in enumerate(givens):
    for c, title in enumerate(titles):
        letter = ord("a") + g*len(titles) + c
        axs[g, c, 0, 0].set_title(chr(letter),
                                  fontsize=letterlabel,
                                  pad = 10,
                                  weight="bold",
                                  loc="left")
        if g == 0:
            axs[0, c, 0, 10].set_title(title,
                                       pad=plotsize*9,
                                       fontsize=letterlabel)
        for a in range(0, nr, step):
            axs[g, c, a, 0].set(yticks=[ylim[1]/2.0], yticklabels=[])
            if c == 0:
                axs[g, 0, a, 0].set_yticklabels([alphas[a]],
                                                fontsize=ticklabel)
        for e in range(0, nc, step):
            axs[g, c, -1, e].set(xticks=[xlim[1]/2.0], xticklabels=[])
            if given == givens[-1]:
                axs[-1, c, -1, e].set_xticklabels([f"{logess[e]:.0f}"],
                                                 fontsize=ticklabel)

if movie:
    fig.text(right_x,
             bottom_y*0.5,
             "t\n0",
             fontsize=biglabel,
             color="grey",
             ha="right")

# Assign axs objects to variables
# (Line2D)

artists = np.empty_like(axs) 
x = np.arange(64)
dummy_y = np.zeros_like(x)
frames = ts
frame0 = ts[-1]

for folder in folders:
    for g, given in enumerate(givens):
        for c, title in enumerate(titles):
            for a, alpha in enumerate(alphas):
                for e, loges in enumerate(logess):
                    ax = axs[g, c, a, e] 
                    artists[g, c, a, e], = ax.plot(x, dummy_y, c="white")

    # Add data and save figure

    dffrqs = np.empty(len(givens), dtype=object) 
    dfmeans = np.empty(len(givens), dtype=object) 
    for g, given in enumerate(givens):
        filelist = glob(os.path.join(folder, given, "*.frq"))
        dffrqs[g] = my.read_files(filelist, movie)
        filelist = glob(os.path.join(folder, given, "*.csv"))
        dfmeans[g] = my.read_files(filelist, movie)

    if movie:
        ani = FuncAnimation(fig,
                            update,
                            frames=frames,
                            fargs=(artists,),
                            blit=True)
        ani.save(f"{file_name}_{folder}.mp4", writer="ffmpeg", fps=10)
    else:
        update(frame0, artists,)
        plt.savefig(f"{file_name}_{folder}.png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
