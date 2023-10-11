#! /usr/bin/env python

from glob import glob
import os
import time

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
this_file = os.path.basename(__file__)
file_name = this_file.split(".")[0]

# Options

factor = "TSselected"
trait = "MimicGrainmean"

movie = False
plotsize = 12

# Add data to figure

def update(t, artists):
    df_t = df[df.Time == t]
    rhos = 1.0 - 1.0/pow(2, df_t.logES)
    T = my.fitness(df_t.a2high, df_t.a2low, df_t.Given, df_t.alpha, rhos)
    R = my.fitness(df_t.a2high, df_t.a2high, df_t.Given, df_t.alpha, rhos)
    P = my.fitness(df_t.a2low, df_t.a2low, df_t.Given, df_t.alpha, rhos)
    S = my.fitness(df_t.a2low, df_t.a2high, df_t.Given, df_t.alpha, rhos)
    e = np.array(R - P)
    e = e.ravel()
    f = np.array(P - S)
    f = f.ravel()
    y0 = np.array(df_t[trait])
    y0 = y0.ravel()
    x0 = np.array(R - (T + S)/2.0)
    x0 = x0.ravel()
    mask = (e > 0.0) & (e < 1.0) & (f > 0.0) & (f < 0.3)
    y = y0[mask]
    x = x0[mask]
    new_offsets = np.dstack((x, y)).reshape(-1, 2)
    artists.set_offsets(new_offsets)

    if movie:
        fig.texts[1].set_text(t)
        return artists, fig.texts[1]
    else:
        return artists

# Data

filelist = glob("*.csv")
if filelist == []:
    print("No *.csv")
    exit()
df = my.read_files(filelist, movie)

ts = df.Time.unique()

# Figure properties

width = plotsize
height = plotsize
biglabel = plotsize*4
ticklabel = plotsize*3
step = 0.5
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Create figure
# Create a figure with a density function of "trait"

fig, ax = plt.subplots(figsize=(width, height))

ax.set_xlim(-my.wmax, my.wmax)
ax.set_ylim(0, my.aBmax)
ax.set_xlabel(factor, fontsize=biglabel)
ax.set_ylabel(trait, fontsize=biglabel)

if movie:
    fig.text(right_x,
             bottom_y*0.5,
             "t\n0",
             fontsize=biglabel,
             color="grey",
             ha="right")

# Assign axs objects to variables
# (Axes)

dummy_x = [-1.0, 0.0]
dummy_y = [0.0, 0.0] 
frames = ts

artists = ax.scatter(dummy_x, dummy_y, s=0.4, c="blue")

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(artists,),
                        blit=True)
    ani.save(f"{factor}_{trait}.mp4", writer="ffmpeg", fps=10)
else:
    update(frames[-1], artists,)
    plt.savefig(f"{factor}_{trait}.png", transparent=False)

plt.close()

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
