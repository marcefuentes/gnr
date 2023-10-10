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
rows = ["pr", "p", "r", "none"]
given = "given100"

movie = False
plotsize = 8
bins = 64

# Add data to figure

def update(t, artists):
    for r, row in enumerate(rows):
        df = dffrqs[r]
        if movie:
            m = df.Time == t
            df = df.loc[m]
        for c, trait in enumerate(traits):
            Z = my.getZ(t, dfmeans[r], trait + "mean")
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
                    artists[r, c, a, e].set_ydata(y)
                    bgcolor = cm.viridis(Z[a, e]/vmaxs[c])
                    artists[r, c, a, e].axes.set_facecolor(bgcolor)
    if movie:
        fig.texts[2].set_text(f"t\n{t}")
    return artists.flatten()

# Data without partner choice or reciprocity

filelist = glob(os.path.join("none", "given000", "*.csv"))
dfsocial = my.read_files(filelist, movie)

filelist = glob(os.path.join("none", given, "*.csv"))
df = my.read_files(filelist, movie)

ts = df.Time.unique()
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(rows)
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

axs = np.empty((len(rows),
                len(titles),
                nr,
                nc),
                dtype=object)

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=len(rows),
                             ncols=len(titles))
for r, row in enumerate(rows):
    for c, title in enumerate(titles):
        grid = outergrid[r, c].subgridspec(nrows=nr,
                                           ncols=nc,
                                           wspace=0,
                                           hspace=0)
        axs[r, c] = grid.subplots()

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

for r, row in enumerate(rows):
    for c, title in enumerate(titles):
        letter = ord("a") + r*len(titles) + c
        axs[r, c, 0, 0].set_title(chr(letter),
                                  fontsize=letterlabel,
                                  pad = 10,
                                  weight="bold",
                                  loc="left")
        if r == 0:
            axs[0, c, 0, 10].set_title(title,
                                       pad=plotsize*9,
                                       fontsize=letterlabel)
        for a in range(0, nr, step):
            axs[r, c, a, 0].set(yticks=[ylim[1]/2.0], yticklabels=[])
            if c == 0:
                axs[r, 0, a, 0].set_yticklabels([alphas[a]],
                                                fontsize=ticklabel)
        for e in range(0, nc, step):
            axs[r, c, -1, e].set(xticks=[xlim[1]/2.0], xticklabels=[])
            if row == rows[-1]:
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

for r, row in enumerate(rows):
    for c, title in enumerate(titles):
        for a, alpha in enumerate(alphas):
            for e, loges in enumerate(logess):
                ax = axs[r, c, a, e] 
                artists[r, c, a, e], = ax.plot(x, dummy_y, c="white")

# Add data and save figure

dffrqs = np.empty(len(rows), dtype=object) 
dfmeans = np.empty(len(rows), dtype=object) 
for r, row in enumerate(rows):
    filelist = glob(os.path.join(row, given, "*.frq"))
    dffrqs[r] = my.read_files(filelist, movie)
    filelist = glob(os.path.join(row, given, "*.csv"))
    dfmeans[r] = my.read_files(filelist, movie)

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(artists,),
                        blit=True)
    ani.save(f"{file_name}_{given}.mp4", writer="ffmpeg", fps=10)
else:
    update(frames[-1], artists,)
    plt.savefig(f"{file_name}_{given}.png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
