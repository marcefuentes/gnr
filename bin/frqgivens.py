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

traits = ["ChooseGrain",
          "MimicGrain"]
titles = ["Sensitivity for\nchoosing partner",
          "Sensitivity for\nmimicking partner"]
vmaxs = [my.a2max, my.a2max]

folders = ["given100", "given095", "given050"]
subfolders = ["p", "r"]

movie = False
plotsize = 8
bins = 64

# Add data to figure

def update(t, artists):
    for f, folder in enumerate(folders):
        for c, trait in enumerate(traits):
            df = dffrqs[f, c]
            if movie:
                m = df.Time == t
                df = df.loc[m]
            Z = my.get_Z(t, dfmeans[f, c], trait + "mean")
            if "Grain" in trait:
                Z = 1.0 - Z
            for a, alpha in enumerate(alphas):
                for l, loges in enumerate(logess):
                    m = (df.alpha == alpha) & (df.logES == loges)
                    d = df.loc[m]
                    freq_a = [col for col in d.columns if re.match(fr"^{trait}\d+$", col)]
                    y = d.loc[:, freq_a]
                    y = y.values[0]
                    y = y.flatten()
                    artists[f, c, a, l].set_ydata(y)
                    bgcolor = cm.viridis(Z[a, l]/vmaxs[c])
                    artists[f, c, a, l].axes.set_facecolor(bgcolor)
    if movie:
        fig.texts[2].set_text(f"t\n{t}")
    return artists.flatten()

# Data

dffrqs = np.empty((len(folders), len(subfolders)), dtype=object)
for f, folder in enumerate(folders):
    for c, subfolder in enumerate(subfolders):
        filelist = glob(os.path.join(folder, subfolder, "*.frq"))
        dffrqs[f, c] = my.read_files(filelist, movie)

dfmeans = np.empty((len(folders), len(subfolders)), dtype=object)
for f, folder in enumerate(folders):
    for c, subfolder in enumerate(subfolders):
        filelist = glob(os.path.join(folder, subfolder, "*.csv"))
        dfmeans[f, c] = my.read_files(filelist, movie)

df = dffrqs[0, 0]
ts = df.Time.unique()
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(folders)
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
xlim = [-2, bins + 1]
ylim = [0, 0.25]
step = int(nr/2)
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Create figure

axs = np.empty((len(folders),
                len(titles),
                nr,
                nc),
                dtype=object)

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=len(folders),
                             ncols=len(titles))

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
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=left_x*0.1,
              y=center_y,
              fontsize=biglabels)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    ax.set(xlim=xlim, ylim=ylim)
    for axis in ["top","bottom","left","right"]:
        ax.spines[axis].set_linewidth(0.1)

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        letter = ord("a") + f*len(titles) + c
        axs[f, c, 0, 0].set_title(chr(letter),
                                  fontsize=plotsize*5,
                                  pad = 11,
                                  weight="bold",
                                  loc="left")
        if f == 0:
            axs[0, c, 0, 10].set_title(title,
                                       pad=plotsize*9,
                                       fontsize=plotsize*5)
        for a in range(0, nr, step):
            axs[f, c, a, 0].set(yticks=[ylim[1]/2.0], yticklabels=[])
            if c == 0:
                axs[f, 0, a, 0].set_yticklabels([alphas[a]],
                                                fontsize=ticklabels)
        for l in range(0, nc, step):
            axs[f, c, -1, l].set(xticks=[xlim[1]/2.0], xticklabels=[])
            if folder == folders[-1]:
                axs[-1, c, -1, l].set_xticklabels([f"{logess[l]:.0f}"],
                                                 fontsize=ticklabels)

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
x = np.arange(64)
dummy_y = np.zeros_like(x)
frames = ts

for f, folder in enumerate(folders):
    for c, title in enumerate(titles):
        for a, alpha in enumerate(alphas):
            for l, loges in enumerate(logess):
                ax = axs[f, c, a, l] 
                artists[f, c, a, l], = ax.plot(x, dummy_y, c="white")

# Add data and save figure

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(artists,),
                        blit=True)
    ani.save(file_name + ".mp4", writer="ffmpeg", fps=10)
else:
    update(frames[-1], artists,)
    plt.savefig(file_name + ".png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
