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

factor = "snowdriftTS_PS"
trait = "MimicGrainmean"

plotsize = 8

# Add data to figure

def update(t, artist, text):
    artist.set_ydata([T[t], R[t], P[t], S[t]])
    bgcolor = cm.viridis(traitcolors[t])
    artist.axes.set_facecolor(bgcolor)
    artist.axes.lines[0].set_ydata([(T[t] + S[t])/2, (T[t] + S[t])/2])
    text.set_text(f"{trait}\nalpha {alphas[t]:.1f}\nlogES  {logESs[t]:2.0f}\na2high {a2highs[t]:.1f}\na2low {a2lows[t]:.1f}")

    return artist, text

# Data

filelist = glob("*.csv")
if filelist == []:
    print("No *.csv")
    exit()
df = my.read_files(filelist, 0)

df = df.sort_values(by=[trait])
a2highs = df.a2high.values.ravel()
a2lows = df.a2low.values.ravel()
givens = df.Given.values.ravel()
alphas = df.alpha.values.ravel()
logESs = df.logES.values.ravel()
traitcolors = df[trait].values.ravel()
rhos = 1.0 - 1.0/pow(2, logESs)
T = my.fitness(a2highs, a2lows, givens, alphas, rhos)
R = my.fitness(a2highs, a2highs, givens, alphas, rhos)
P = my.fitness(a2lows, a2lows, givens, alphas, rhos)
S = my.fitness(a2lows, a2highs, givens, alphas, rhos)
mask = (T > R) & (R > P) & (P < S) & (P == 0) & (T + S > 2*R)
T = T[mask]
R = R[mask]
P = P[mask]
S = S[mask]
traitcolors = traitcolors[mask]

# Figure properties

width = plotsize
height = plotsize
biglabel = plotsize*4
ticklabel = plotsize*3
step = 0.5
xlim = [0, 5]
ylim = [0.0, my.wmax]
xticks = [1, 2, 3, 4]
yticks = [0.0, my.wmax/2, my.wmax]
xticklabels = ["T", "R", "P", "S"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Create figure

fig, ax = plt.subplots(figsize=(width, height))

text = ax.text(0.95,
               0.82,
               "dummy",
               fontsize=ticklabel/2,
               color="white",
               ha="right",
               fontname="monospace",
               transform=ax.transAxes)

ax.set(xlim=xlim, ylim=ylim)
ax.set(xticks=xticks, yticks=yticks)
ax.set_xticklabels(xticklabels, fontsize=ticklabel)
ax.set_yticklabels(yticks, fontsize=ticklabel)
ax.axhline(y=my.wmax/2, color="white", linestyle="--", linewidth=1)

# Assign axs objects to variables
# (Line2D)

xaxis = [1, 2, 3, 4]
dummy_y = [0.5, 0.5, 0.5, 0.5]
frames = range(len(traitcolors))

artist, = ax.plot(xaxis,
                  dummy_y,
                  color="white",
                  marker="o",
                  markerfacecolor="white")

ani = FuncAnimation(fig,
                    update,
                    frames=frames,
                    fargs=(artist, text,),
                    blit=True)
ani.save(f"{trait}_trps.mp4", writer="ffmpeg", fps=10)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
