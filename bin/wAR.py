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

numi = 1025 # Number of inner plot values
numo = 11  # Number of outer plot values

movie = False
if movie:
    givens = np.linspace(0.0, 1.0, num=41)
else:
    givens = [0.95]
plotsize = 12

# Add data to figure

def update(given, artists):
    for y, alpha in enumerate(alphas):
        AA = np.full([numi, numi], alpha)
        for x, rho in enumerate(rhos):
            RR = np.full([numi, numi], rho)
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
            Z[XX >= YY] = 0.0
            artists[y, x].set_array(Z)
    if movie:
        fig.texts[2].set_text(f"{given:.2f}")
    return artists.flatten()

# Data

alphas = np.linspace(my.alphamax, my.alphamin, num=numo)
logess = np.linspace(my.logesmin, my.logesmax, num=numo)
rhos = 1.0 - 1.0/pow(2, logess)
aBys = np.linspace(my.aBmax, 0.0, num=numi)
aBxs = np.linspace(0.0, my.aBmax, num=numi)
XX, YY = np.meshgrid(aBxs, aBys)

# Figure properties

width = plotsize
height = plotsize
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
biglabel = plotsize*4
ticklabel = plotsize*3
extent = 0, numi, 0, numi
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
    axs[y, 0].set_ylabel(f"{alphas[y]:.1f}",
                         rotation="horizontal",
                         horizontalalignment="right",
                         verticalalignment="center",
                         fontsize=ticklabel)
for x in range(0, numo, step):
    axs[-1, x].set_xlabel(f"{logess[x]:.0f}",
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
dummy_Z = np.full(XX.shape, 0.0)
frames = givens
frames0 = frames[0]

for y, alpha in enumerate(alphas):
    for x, rho in enumerate(rhos):
        artists[y, x] = axs[y, x].imshow(dummy_Z,
                                         extent=extent,
                                         vmin=0.0,
                                         vmax=my.wmax,
                                         aspect="auto",
                                         interpolation="nearest")

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
