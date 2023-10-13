#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.transforms
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
this_file = os.path.basename(__file__)
file_name = this_file.split(".")[0]

# Options

traits = ["ChooseGrainmean",
          "MimicGrainmean",
          "ImimicGrainmean"]
titles = ["Partner choice",
          "Direct reciprocity",
          "Indirect reciprocity"]

numa2 = 64
plotsize = 6
r = 0.1

# Data

filelist = glob("../../none/given000/*.csv")
df = my.read_files(filelist, False)
t = df.Time.iloc[0]
a2socials = my.getZ(t, df, "a2Seenmean")
wsocials = my.getZ(t, df, "wmean")

folder = os.path.basename(os.getcwd())
filelist = glob(os.path.join("../../none", folder, "*.csv"))
df = my.read_files(filelist, False)
a2privates = my.getZ(t, df, "a2Seenmean")
wprivates = my.getZ(t, df, "wmean")

filelist = glob("*.csv")
df = my.read_files(filelist, False)
alphas = np.sort(pd.unique(df.alpha))[::-1]
logess = np.sort(pd.unique(df.logES))
nr = len(alphas)
nc = len(logess)
rhos = 1. - 1./pow(2., logess)
given = df.Given.iloc[0]
xaxis = np.linspace(0.01, my.a2max - 0.01, num=numa2)
Z = np.zeros((len(traits), len(alphas), len(rhos)))
for c, trait in enumerate(traits):
    Z[c] = my.getZ(t, df, trait)
    if "Grain" in trait:
        Z[c] = 1.0 - Z[c]

# Figure properties

width = plotsize*len(titles)
height = plotsize
xlabel = "Substitutability of $\it{B}$"
ylabel = "Influence of $\it{B}$"
biglabel = plotsize*4
letterlabel = plotsize*3
ticklabel = plotsize*3
xlim = [0., my.a2max]
ylim = [0., 1.]
step = int(nr/2)
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Create figure

fig = plt.figure(figsize=(width, height))
outergrid = fig.add_gridspec(nrows=1,
                             ncols=len(titles),
                             left=0.1,
                             right=0.9,
                             top=0.85,
                             bottom=0.17)

axs = np.empty((len(titles),
                nr,
                nc),
               dtype=object)

for c, title in enumerate(titles):
    grid = outergrid[c].subgridspec(nrows=nr,
                                    ncols=nc,
                                    wspace=0,
                                    hspace=0)
    axs[c] = grid.subplots()

left_x = axs[0, 0, 0].get_position().x0
right_x = axs[-1, -1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0, 0].get_position().y1
bottom_y = axs[-1, -1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.2,
              fontsize=biglabel)
fig.supylabel(ylabel,
              x=left_x*0.4,
              y=center_y,
              fontsize=biglabel)

for ax in fig.get_axes():
    ax.set(xticks=[], yticks=[])
    ax.set(xlim=xlim, ylim=ylim)
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(0.1)

for c, title in enumerate(titles):
    letter = ord("a") + c
    axs[c, 0, 0].set_title(chr(letter),
                           fontsize=letterlabel,
                           pad = 11,
                           weight="bold",
                           loc="left")
    axs[c, 0, 10].set_title(title,
                            pad=plotsize*5,
                            fontsize=biglabel)
    for a in range(0, nr, step):
        axs[c, a, 0].set(yticks=[ylim[1]/2.0], yticklabels=[])
        if c == 0:
            axs[0, a, 0].set_yticklabels([alphas[a]],
                                         fontsize=ticklabel)
    for e in range(0, nc, step):
        axs[c, -1, e].set(xticks=[xlim[1]/2.0], xticklabels=[])
        axs[c, -1, e].set_xticklabels([f"{logess[e]:.0f}"],
                                      fontsize=ticklabel)

# Add data

for a, alpha in enumerate(alphas):
    for e, rho in enumerate(rhos):

        a2s = np.full(xaxis.shape, a2privates[a, e])

        ax = axs[0, a, e]
        ii = my.fitness(xaxis, xaxis, given, alpha, rho)
        ji = my.fitness(a2s, xaxis, given, alpha, rho)
        jj = my.fitness(a2s, a2s, given, alpha, rho)
        ij = my.fitness(xaxis, a2s, given, alpha, rho)
        y = (jj - ii*r - ji*(1.0 - r))/((1.0 - r)*(ii - ij - ji + jj))
        mask = xaxis < a2privates[a, e]
        y[mask] = np.nan
        ax.plot(xaxis, y, color="white", linewidth=0.7)
        #ax.plot(a2privates[a, e],
        #        wprivates[a, e],
        #        marker="o",
        #        markersize=1.5,
        #        color="white")

        ax = axs[1, a, e]
        y = (my.repeats*my.cost + jj - ji)/(my.repeats*ii - ij - ji + 2.*jj - my.repeats*jj) 
        ax.plot(xaxis, y, color="white", linewidth=0.7)
        ax.plot(a2privates[a, e],
                wprivates[a, e],
                marker="o",
                markersize=1.5,
                color="white")

        ax = axs[2, a, e]
        ax.plot(xaxis, y, color="white", linewidth=0.7)
        ax.plot(a2privates[a, e],
                wprivates[a, e],
                marker="o",
                markersize=1.5,
                color="white")

        for c, title in enumerate(titles):
            color = cm.viridis(Z[c, a, e]/my.a2max)
            axs[c, a, e].set_facecolor(color)

# Finish

plt.savefig(file_name + ".png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
