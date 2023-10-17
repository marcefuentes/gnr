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

numi = 257 # Number of inner plot values
numo = 11  # Number of outer plot values

movie = False
if movie:
    givens = np.linspace(0.0, 1.0, num=41)
else:
    givens = [1.00]
plotsize = 12

# Add data to figure

def update(given, artists):
    for y, a2y in enumerate(a2ys):
        for x, a2x in enumerate(a2xs[:numo-y-1]):
            XX = np.full(AA.shape, a2x)
            YY = np.full(AA.shape, a2y)
            T = my.fitness(YY, XX, given, AA, RR)
            R = my.fitness(YY, YY, given, AA, RR)
            P = my.fitness(XX, XX, given, AA, RR)
            S = my.fitness(XX, YY, given, AA, RR)
            w = my.eqw(T, R, P, S)
            T = my.fitness(YY, XX, 0.0, AA, RR)
            R = my.fitness(YY, YY, 0.0, AA, RR)
            P = my.fitness(XX, XX, 0.0, AA, RR)
            S = my.fitness(XX, YY, 0.0, AA, RR)
            wsocial = my.eqw(T, R, P, S)
            Z = wsocial - w
            artists[y, x].set_array(Z)
    if movie:
        fig.texts[2].set_text(f"{given:.2f}")
    return artists.flatten()

# Data

alphas = np.linspace(my.alphamax, my.alphamin, num=numi)
logess = np.linspace(my.logesmin, my.logesmax, num=numi)
rhos = 1.0 - 1.0/pow(2, logess)
RR, AA = np.meshgrid(rhos, alphas)
a2xs = np.linspace(0.0, my.a2max, num=numo)
a2ys = np.linspace(my.a2max, 0.0, num=numo)

# Figure properties

width = plotsize
height = plotsize
xlabel = "Effort to get $\it{B}$"
ylabel = "Effort to get $\it{B}$"
biglabel = plotsize*4
ticklabel = plotsize*3
step = int(numo/2)
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Create figure

fig = plt.figure(figsize=(width, height))
grid = fig.add_gridspec(nrows=numo,
                        ncols=numo,
                        left=0.22,
                        right=0.9,
                        top=0.86,
                        bottom=0.176,
                        wspace=0,
                        hspace=0)

axs = grid.subplots()

left_x = axs[0, 0].get_position().x0
right_x = axs[-1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0].get_position().y1
bottom_y = axs[-1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.2,
              fontsize=biglabel)
fig.supylabel(ylabel,
              x=left_x*0.2,
              y=center_y,
              fontsize=biglabel)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    for axis in ["top","bottom","left","right"]:
        ax.spines[axis].set_linewidth(0.2)
for y in range(0, numo, step):
    axs[y, 0].set_ylabel(f"{a2ys[y]:.1f}",
                         rotation="horizontal",
                         horizontalalignment="right",
                         verticalalignment="center",
                         fontsize=ticklabel)
for x in range(0, numo, step):
    axs[-1, x].set_xlabel(f"{a2xs[x]:.1f}",
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
dummy_Z = np.full((numi, numi, 4), my.colormap["greyTS"])
frames = givens
frames0 = frames[0]

for y, a2y in enumerate(a2ys):
    for x, a2x in enumerate(a2xs):
        artists[y, x] = axs[y, x].imshow(dummy_Z,
                                         vmin=0.0,
                                         vmax=my.wmax,
                                         aspect="auto")

# Add data and save figure

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(artists,),
                        blit=True)
    ani.save(f"{file_name}.mp4", writer="ffmpeg", fps=10)
else:
    update(frames0, artists,)
    plt.savefig(f"{file_name}.png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
