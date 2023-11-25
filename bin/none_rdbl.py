#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import matplotlib.transforms
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
this_file = os.path.basename(__file__)
file_name = this_file.split(".")[0]

# Options

titles = ["Production of $\it{B}$",
          "Byproduct help",
          "Fitness",
          "Fitness deficit"]
givens = ["given100", "given095", "given050", "given000"]
movie = False
plotsize = 4

# Add data to figure

def update(t, artists):
    Zsocial = my.getZ(t, dfs[-1], "wmean")/my.wmax 
    for g, given in enumerate(givens):
        Z = my.getZ(t, dfs[g], "a2Seenmean")
        artists[g, 0].set_array(Z)
        given_float = float(given[-3:])/100
        Z = Z*given_float
        artists[g, 1].set_array(Z)
        Z = my.getZ(t, dfs[g], "wmean")/my.wmax
        artists[g, 2].set_array(Z)
        Z = Z - Zsocial
        artists[g, 3].set_array(Z)
    return artists.flatten()

# Data

dfs = np.empty(len(givens), dtype=object)
for g, given in enumerate(givens):
    filelist = glob(os.path.join("none", given, "*.csv"))
    dfs[g] = my.read_files(filelist, movie)

df = dfs[0]
ts = df.Time.unique()
nr = df.alpha.nunique()
nc = df.logES.nunique()
rows = givens

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(rows)
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
biglabel = plotsize*7
letterlabel = plotsize*6
ticklabel = plotsize*5
xticks = [0, nc/2 - 0.5, nc - 1]
yticks = [0, nr/2 - 0.5, nr - 1]
xmin = df.logES.min()
xmax = df.logES.max()
ymin = df.alpha.min()
ymax = df.alpha.max()
xticklabels = [f"{xmin:.0f}",
               f"{(xmin + xmax)/2.:.0f}",
               f"{xmax:.0f}"]
yticklabels = [f"{ymax:.1f}",
               f"{(ymin + ymax)/2.:.1f}",
               f"{ymin:.1f}"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Create figure

fig, axs = plt.subplots(nrows=len(rows),
                        ncols=len(titles),
                        figsize=(width, height))

left_x = axs[0, 0].get_position().x0
right_x = axs[-1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0].get_position().y1
bottom_y = axs[-1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y - 1.2/height,
              fontsize=biglabel)
fig.supylabel(ylabel,
              x=left_x - 1.45/width,
              y=center_y,
              fontsize=biglabel)

letterposition = 1.035
for i, ax in enumerate(fig.get_axes()):
    ax.set(xticks=xticks, yticks=yticks)
    ax.set(xticklabels=[], yticklabels=[])
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(0.1)
    letter = ord("a") + i
    ax.text(0,
            letterposition,
            chr(letter),
            transform=ax.transAxes,
            fontsize=letterlabel,
            weight="bold")
for r, row in enumerate(rows):
    axs[r, 0].set_yticklabels(yticklabels, fontsize=ticklabel)
for c, title in enumerate(titles):
    axs[0, c].set_title(title, pad=plotsize*10, fontsize=letterlabel)
    axs[-1, c].set_xticklabels(xticklabels,
                               fontsize=ticklabel)
if movie:
    fig.text(right_x,
             bottom_y*0.5,
             "t\n0",
             fontsize=biglabel,
             color="grey",
             ha="right")

# Assign axs objects to variables
# (AxesImage)

artists = np.empty_like(axs) 
dummy_Z = np.zeros((nr, nc))
frames = ts

for r, row in enumerate(rows):
    for c, title in enumerate(titles):
        artists[r, c] = axs[r, c].imshow(dummy_Z,
                                         cmap="RdBu_r",
                                         vmin=-1,
                                         vmax=1)

axins = inset_axes(axs[-1, -1],
                   width="5%",
                   height="100%",
                   loc="upper right",
                   bbox_to_anchor=(1220, 540, 300, 500),
                   borderpad=0)
cbar = fig.colorbar(artists[-1, -1],
                    cax=axins,
                    ticks=[-1, 0, 1])
cbar.ax.tick_params(labelsize=ticklabel)
cbar.outline.set_linewidth(0.1)

# Save figure

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(artists,),
                        blit=True)
    ani.save(f"{file_name}.mp4", writer="ffmpeg", fps=10)
else:
    update(frames[-1], artists,)
    plt.savefig(f"{file_name}.png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
