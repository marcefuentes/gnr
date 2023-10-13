#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib.animation import FuncAnimation
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
vmaxs = [my.a2max,
         my.a2max,
         my.wmax,
         my.wmax]

givens = np.linspace(0., 1., num=42)
movie = True
plotsize = 4

# Add data to figure

def update(given, artists):
    a2 = my.a2eq(given, AA, RR)
    w = my.fitness(a2, a2, given, AA, RR)
    dif = wsocial - w
    artists[0].set_array(a2)
    artists[1].set_array(a2*given)
    artists[2].set_array(w)
    artists[3].set_array(dif)
    if movie:
        fig.texts[2].set_text(f"Given {given*100:.0f}%")
    return artists.flatten()

# Data

nr = 513
nc = nr
alphas = np.linspace(my.alphamax, my.alphamin, num=nr)
logess = np.linspace(my.logesmin, my.logesmax, num=nc)
rhos = 1. - 1./pow(2., logess)
RR, AA = np.meshgrid(rhos, alphas)
a2social = my.a2eq(0., AA, RR)
wsocial = my.fitness(a2social, a2social, 0., AA, RR)

# Figure properties

width = plotsize*len(titles) - 2.
height = plotsize
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
biglabel = plotsize*6
letterlabel = plotsize*5
ticklabel = plotsize*4
xticks = [0, nc/2 - 0.5, nc - 1]
yticks = [0, nr/2 - 0.5, nr - 1]
xmin = my.logesmin
xmax = my.logesmax
ymin = my.alphamin
ymax = my.alphamax
xticklabels = [f"{xmin:.0f}",
               f"{(xmin + xmax)/2.:.0f}",
               f"{xmax:.0f}"]
yticklabels = [f"{ymax:.1f}",
               f"{(ymin + ymax)/2.:.1f}",
               f"{ymin:.1f}"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Create figure

fig, axs = plt.subplots(nrows=1,
                        ncols=len(titles),
                        figsize=(width, height))
plt.subplots_adjust(top=0.80, bottom=0.25, wspace=0.3)

left_x = axs[0].get_position().x0
right_x = axs[-1].get_position().x1
center_x = (left_x + right_x) / 2.
top_y = axs[0].get_position().y1
bottom_y = axs[0].get_position().y0
center_y = (top_y + bottom_y) / 2.
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.2,
              fontsize=biglabel)
fig.supylabel(ylabel,
              x=left_x*0.4,
              y=center_y,
              fontsize=biglabel)

ox = 0/72.
oy = 0/72.
offset = matplotlib.transforms.ScaledTranslation(ox, oy, fig.dpi_scale_trans)

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
axs[0].set_yticklabels(yticklabels, fontsize=ticklabel)
for c, title in enumerate(titles):
    axs[c].set_title(title, pad=plotsize*7, fontsize=letterlabel)
    axs[c].set_xticklabels(xticklabels,
                               fontsize=ticklabel)
    for label in axs[c].xaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)

if movie:
    fig.text(right_x,
             bottom_y*0.2,
             f"Given",
             fontsize=biglabel,
             color="grey",
             ha="right")

# Assign axs objects to variables
# (AxesImage)

artists = np.empty_like(axs) 
dummy_Z = np.zeros((nr, nc))
frames = givens

for c, title in enumerate(titles):
    artists[c] = axs[c].imshow(dummy_Z,
                               vmin=0,
                               vmax=vmaxs[c])

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
