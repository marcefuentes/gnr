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
title = "Reciprocity"
vmax = my.aBmax

movie = False
plotsize = 12

# Add data to figure

def init(artists):

    for a, alpha in enumerate(alphas):
        AA = np.full((numi, numi), alpha)
        for r, rho in enumerate(rhos):
            RR = np.full((numi, numi), rho)
            T = my.fitness(YY, XX, given, AA, RR)
            R = my.fitness(YY, YY, given, AA, RR)
            P = my.fitness(XX, XX, given, AA, RR)
            S = my.fitness(XX, YY, given, AA, RR)
            #Ma = np.maximum.reduce([T, R, P, S])
            #Mi = np.minimum.reduce([T, R, P, S])
            #Tn = (T - Mi)/(Ma - Mi)
            #Rn = (R - Mi)/(Ma - Mi)
            #Pn = (P - Mi)/(Ma - Mi)
            #Sn = (S - Mi)/(Ma - Mi)
            Tn, Rn, Pn, Sn = T, R, P, S
            y = np.stack((Tn, Rn, Pn, Sn), axis=-1)
            linecolor = np.full(highs.shape, "white")
            red = np.full(highs.shape, "red")
            m = lows > highs
            linecolor[m] = red[m]

            for (high, low, i), _ in np.ndenumerate(y):
                artists[a, r, high, low].set_ydata(y[high, low])
                lcolor = linecolor[high, low] 
                artists[a, r, high, low].set_color(lcolor)
                artists[a, r, high, low].set_markerfacecolor(lcolor)

    return artists.flatten()

def update(t, artists):
    for a, alpha in enumerate(alphas):
        for r, loges in enumerate(logess):
            Z = my.getZd(t, df, alpha, loges, trait)
            #a2lows = my.getZd(t, df, alpha, loges, "a2low")
            #Z = Z - a2lows
            if "Grain" in trait:
                Z = 1.0 - Z
            if "gain" in title:
                wnull = my.getZ(t, df, "wmean")
                Z = Z - wnull
            if "deficit" in title:
                wsocial = my.getZd(t, dfsocial, alpha, loges, trait)
                wnone = my.getZd(t, dfnone, alpha, loges, trait)
                Z = wsocial - wnone
            for y, aBy in enumerate(aBys):
                for x, aBx in enumerate(aBxs):
                    if aBy <= aBx:
                        artists[a, r, y, x].axes.set_facecolor("white")
                    else:
                        bgcolor = cm.viridis(Z[y, x]/vmax)
                        artists[a, r, y, x].axes.set_facecolor(bgcolor)
    if movie:
        fig.texts[2].set_text(t)
    return artists.flatten()

# Data

if "deficit" in title:
    filelist = glob("../../none/given000/*.csv")
    if filelist == []:
        print("No *.csv")
        exit()
    dfsocial = my.read_files(filelist, movie)
    filelist = glob("../../none/given100/*.csv")
    if filelist == []:
        print("No *.csv")
        exit()
    dfnone = my.read_files(filelist, movie)

filelist = glob("*.csv")
if filelist == []:
    print("No *.csv")
    exit()
df = my.read_files(filelist, movie)

ts = df.Time.unique()
alphas = np.sort(df.alpha.unique())[::-1]
logess = np.sort(df.logES.unique())
rhos = 1.0 - 1.0/pow(2, logess)
given = df.Given.iloc[0]
numo = len(alphas)
numi = df.a2low.nunique() + 1
aBxs = np.linspace(0.0, my.aBmax, num=numi)
aBys = np.linspace(my.aBmax, 0.0, num=numi)
XX, YY = np.meshgrid(aBxs, aBys)
numi2 = int((numi + 1)/2)

# Figure properties

width = plotsize
height = plotsize
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
biglabel = plotsize*4
ticklabel = plotsize*3
step = int(numo/2)
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Create figure

axs = np.empty([numo,
                numo,
                numi,
                numi],
               dtype=object)

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=numo,
                             ncols=numo,
                             left=0.22,
                             right=0.9,
                             top=0.86,
                             bottom=0.176,
                             wspace=0,
                             hspace=0)

for a, alpha in enumerate(alphas):
    for r, rho in enumerate(rhos):
        grid = outergrid[a, r].subgridspec(nrows=numi,
                                           ncols=numi,
                                           wspace=0,
                                           hspace=0)
        axs[a, r] = grid.subplots()

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
    for axis in ["top","bottom","left","right"]:
        ax.spines[axis].set_linewidth(0.1)
for y in range(0, numo, step):
    axs[y, 0, numi2, 0].set_ylabel(f"{alphas[y]:.1f}",
                                   rotation="horizontal",
                                   horizontalalignment="right",
                                   verticalalignment="center",
                                   fontsize=ticklabel)
for x in range(0, numo, step):
    axs[-1, x, -1, numi2].set_xlabel(f"{logess[x]:.0f}",
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
dummy_y = np.zeros_like(xaxis)
frames = ts
frames0 = frames[0]

for a, alpha in enumerate(alphas):
    for r, rho in enumerate(rhos):
        for y in range(numi):
            for x in range(numi):
                ax = axs[a, r, y, x] 
                artists[a, r, y, x], = ax.plot(xaxis,
                                               dummy_y,
                                               linewidth=0,
                                               marker="o",
                                               markersize=plotsize/300)

# Add data and save figure

#init(artists,)

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
