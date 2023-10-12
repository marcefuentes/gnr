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

numi = 11 # Number of inner plot values
numo = 11  # Number of outer plot values
movie = False
plotsize = 48

# Add data to figure

def update(given, artists):
    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            for k, y in enumerate(ys):
                for l, x in enumerate(xs):
                    if y > x:
                        T = my.fitness(y, x, given, alpha, rho)
                        R = my.fitness(y, y, given, alpha, rho)
                        P = my.fitness(x, x, given, alpha, rho)
                        S = my.fitness(x, y, given, alpha, rho)
                        Z = my.gamecolors(T, R, P, S)
                        artists[i, j, k, l].set_ydata([T, R, P, S])
                        artists[i, j, k, l].axes.set_facecolor(Z)
    if movie:
        fig.texts[3].set_text(t)
    return artists.flatten()

# Data

if movie:
    givens = np.linspace(0.0, 1.0, num=41)
else:
    g = float(os.getcwd()[-3:])/100
    print(f"Given: {g}")
alphas = np.linspace(my.alphamax, my.alphamin, num=numo)
logess = np.linspace(my.logesmin, my.logesmax, num=numo)
rhos = 1.0 - 1.0/pow(2, logess)
ys = np.linspace(my.aBmax, 0.0, num=numi)
xs = np.linspace(0.0, my.aBmax, num=numi)
numi2 = int(numi/2)

# Figure properties

width = plotsize
height = plotsize
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
biglabel = plotsize*4
ticklabel = plotsize*3
step = int(numo/2)
xlim = [0, 5]
ylim = [-0.1, my.wmax]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Create figure

axs = np.empty((numo,
                numo,
                numi,
                numi),
               dtype=object)

fig = plt.figure(figsize=(width, height))
grid = fig.add_gridspec(nrows=numo,
                        ncols=numo,
                        left=0.22,
                        right=0.9,
                        top=0.86,
                        bottom=0.176,
                        wspace=0,
                        hspace=0)

for i in range(numo):
    for j in range(numo):
        innergrid = grid[i, j].subgridspec(nrows=numi,
                                           ncols=numi,
                                           wspace=0,
                                           hspace=0)
        axs[i, j] = innergrid.subplots()

left_x = axs[0, 0, 0, 0].get_position().x0
right_x = axs[-1, -1, -1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0, 0, 0].get_position().y1
bottom_y = axs[-1, -1, -1, -1].get_position().y0
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
    ax.set(xlim=xlim, ylim=ylim)
    for axis in ["top","bottom","left","right"]:
        ax.spines[axis].set_linewidth(0.1)
for i in range(0, numo, step):
    axs[i, 0, numi2, 0].set_ylabel(f"{alphas[i]:.1f}",
                                   rotation="horizontal",
                                   horizontalalignment="right",
                                   verticalalignment="center",
                                   fontsize=ticklabel)
for j in range(0, numo, step):
    axs[-1, j, -1, numi2].set_xlabel(f"{logess[j]:.0f}",
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
xaxis = [1, 2, 3, 4]
dummy_y = [0.5, 0.5, 0.5, 0.5]
if movie:
    frames = givens

for i in range(numo):
    for j in range(numo):
        for k in range(numi):
            for l in range(numi):
                ax = axs[i, j, k, l] 
                artists[i, j, k, l], = ax.plot(xaxis,
                                               dummy_y,
                                               linewidth=0.3,
                                               color="white",
                                               marker="o",
                                               markerfacecolor="white",
                                               markersize=plotsize/40)

# Add data and save figure

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(artists,),
                        blit=True)
    ani.save(f"{file_name}.mp4", writer="ffmpeg", fps=10)
else:
    update(g, artists,)
    plt.savefig(f"{file_name}.png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
