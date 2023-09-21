#! /usr/bin/env python

import os
import time

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

import mymodule as my

start_time = time.perf_counter()
this_file = os.path.basename(__file__)
file_name = this_file.split(".")[0]

# Options

aBx = 0.20
aBy = 0.40
numi = 1025 # Number of inner plot values

movie = False
if movie:
    givens = np.linspace(0.0, 1.0, num=21)
else:
    givens = [0.95]
plotsize = 8

# Add data to figure

def update(given, artists):
    T = my.fitness(YY, XX, given, AA, RR)
    R = my.fitness(YY, YY, given, AA, RR)
    P = my.fitness(XX, XX, given, AA, RR)
    S = my.fitness(XX, YY, given, AA, RR)
    Z = my.gamecolors(T, R, P, S)
    artist.set_array(Z)
    if movie:
        fig.texts[0].set_text(f"{given:.2f}")
    return artist,

# Data

xmin = my.logesmin
xmax = my.logesmax
ymin = my.alphamin
ymax = my.alphamax
alphas = np.linspace(ymax, ymin, num=numi)
logess = np.linspace(xmin, xmax, num=numi)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
XX = np.full(AA.shape, aBx)
YY = np.full(AA.shape, aBy)

# Figure properties

xticks = [0, numi/2.0, numi-0.5]
yticks = [0, numi/2.0, numi-0.5]
xticklabels = [f"{round(xmin, 3):.0f}",
               f"{round((xmax + xmin)/2.0, 3):.0f}",
               f"{round(xmax, 3):.0f}"]
yticklabels = [f"{round(ymax, 3):.1f}",
               f"{round((ymax + ymin)/2.0, 3):.1f}",
               f"{round(ymin, 3):.1f}"]
width = plotsize
height = plotsize
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
biglabel = plotsize*4
ticklabel = plotsize*3
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Create figure

fig, ax = plt.subplots(figsize=(width, height))
plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

ax.set_xlabel(xlabel,
              labelpad=10,
              fontsize=biglabel)
ax.set_ylabel(ylabel,
              labelpad=10,
              fontsize=biglabel)
ax.set(xticks=xticks, yticks=yticks)
ax.set_xticklabels(xticklabels, fontsize=ticklabel)
ax.set_yticklabels(yticklabels, fontsize=ticklabel)

if movie:
    fig.text(0.9,
             0.92,
             "t\n0",
             fontsize=biglabel,
             color="grey",
             ha="right")

# Assign axs objects to variables
# (AxesImage)

dummy_Z = np.full((numi, numi, 4), my.colormap["greyTS"])
frames = givens
frames0 = frames[0]

artist = ax.imshow(dummy_Z)

# Add data and save figure

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(artist,),
                        blit=True)
    ani.save(f"{file_name}.mp4", writer="ffmpeg", fps=10)
else:
    update(frames0, artist,)
    plt.savefig(f"{file_name}.png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
