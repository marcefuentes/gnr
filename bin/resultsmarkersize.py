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

trait = "MimicGrainmean"
title = trait
if "wmean" in trait:
    vmax = my.wmax
else:
    vmax = my.aBmax

movie = False
plotsize = 48

# Add data to figure

def init(artists):
    for i, alpha in enumerate(alphas):
        for j, rho in enumerate(rhos):
            for k, y in enumerate(ys):
                for l, x in enumerate(xs):
                    if y > x:
                        R0 = R[i, j, k, l]
                        P0 = P[i, j, k, l]
                        T0 = T[i, j, k, l]
                        S0 = S[i, j, k, l]
                        condition = (T0 > R0) & (R0 > P0) & (P0 < S0) & (T0 + S0 > R0*2.0)
                        if condition:
                            size = [(S0 - P0)*40]
                            artists[i, j, k, l].set_sizes(size)
    return artists.flatten()

def update(t, artists):
    for i, alpha in enumerate(alphas):
        for j, loges in enumerate(logess):
            Z = my.getZd(t, df, alpha, loges, trait)
            if "Grain" in trait:
                Z = 1.0 - Z
            if "gain" in title:
                wnull = my.getZ(t, df, "wmean")
                Z = Z - wnull
            if "deficit" in title:
                Zsocial = my.getZd(t, dfsocial, alpha, loges, trait)
                Z = Zsocial - Z
            if "a2" in trait:
                Z = (Z - a2lows)/(a2highs - a2lows)
            for k, y in enumerate(ys):
                for l, x in enumerate(xs):
                    if y > x:
                        R0 = R[i, j, k, l]
                        P0 = P[i, j, k, l]
                        T0 = T[i, j, k, l]
                        S0 = S[i, j, k, l]
                        condition = (T0 > R0) & (R0 > P0) & (P0 < S0) & (T0 + S0 > R0*2.0)
                        if condition:
                            bgcolor = cm.viridis(Z[k, l]/vmax)
                            artists[i, j, k, l].axes.set_facecolor(bgcolor)
    if movie:
        fig.texts[3].set_text(t)
    return artists.flatten()

# Data

if "deficit" in title:
    filelist = glob("../../none/given000/*.csv")
    if filelist == []:
        print("No *.csv")
        exit()
    dfsocial = my.read_files(filelist, movie)

filelist = glob("*.csv")
if filelist == []:
    print("No *.csv")
    exit()
df = my.read_files(filelist, movie)

ts = df.Time.unique()
given = df.Given.iloc[0]
alphas = np.sort(df.alpha.unique())[::-1]
logess = np.sort(df.logES.unique())
rhos = 1.0 - 1.0/pow(2, logess)
numo = len(alphas)
if "a2" in trait:
    a2highs = my.getZd(ts[0], df, alphas[0], logess[0], "a2high")
    a2lows = my.getZd(ts[0], df, alphas[0], logess[0], "a2low")
ys = np.sort(df.a2high.unique())[::-1]
xs = np.sort(df.a2low.unique())
ysplus = np.append(ys, 0.0)
xsplus = np.append(xs, 1.0)
XX, YY = np.meshgrid(xsplus, ysplus)
numi = len(ys) + 1
numi2 = int(numi/2)

T = np.zeros((numo, numo, numi, numi))
R = np.zeros((numo, numo, numi, numi))
P = np.zeros((numo, numo, numi, numi))
S = np.zeros((numo, numo, numi, numi))
for i, alpha in enumerate(alphas):
    for j, rho in enumerate(rhos):
        T[i, j] = my.fitness(YY, XX, given, alpha, rho)
        R[i, j] = my.fitness(YY, YY, given, alpha, rho)
        P[i, j] = my.fitness(XX, XX, given, alpha, rho)
        S[i, j] = my.fitness(XX, YY, given, alpha, rho)

# Figure properties

width = plotsize
height = plotsize
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
biglabel = plotsize*4
ticklabel = plotsize*3
step = int(numo/2)
xlim = [0, 1]
ylim = [0, 1]
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
fig.text(center_x,
         top_y*1.05,
         title,
         fontsize=biglabel,
         ha="center")

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
# (Axes)

artists = np.empty_like(axs) 
xaxis = [0.5]
dummy_y = [0.5]
frames = ts

for i in range(numo):
    for j in range(numo):
        for k in range(numi):
            for l in range(numi):
                ax = axs[i, j, k, l] 
                artists[i, j, k, l] = ax.scatter(xaxis,
                                                 dummy_y,
                                                 color="white",
                                                 marker="o")

# Add data and save figure

init(artists,)

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(artists,),
                        blit=True)
    ani.save(f"{title}_size.mp4", writer="ffmpeg", fps=10)
else:
    update(frames[-1], artists,)
    plt.savefig(f"{title}_size.png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
