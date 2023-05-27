#! /usr/bin/env python

import os
import time

from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

# Options

titles = ['Games',
          '$\it{R}$ - $\it{P}$',
          '$\it{T}$ + $\it{S}$ - 2$\it{R}$']
givens = [1.0, 0.95, 0.50]
distance = 0.01
nr = 21
nc = nr
nframes = 21
movie = True
plotsize = 4

# Add data to figure

def update(low, artists):
    for g, given in enumerate(givens):
        lows = np.full(AA.shape, low)
        highs = lows + distance
        T = my.fitness(highs, lows, given, AA, RR)
        R = my.fitness(highs, highs, given, AA, RR)
        P = my.fitness(lows, lows, given, AA, RR)
        S = my.fitness(lows, highs, given, AA, RR)

        Z = my.gamecolors(T, R, P, S)
        artists[g, 0].set_array(Z)
        N = my.nodilemmacolors(T, R, P, S)

        Z = 2.0*R - P - T 
        # m = R < P
        #Z[m] = T[m] - R[m]
        #G = np.full([*Z.shape, 4], my.colormap['transparent'])
        artists[g, 1].set_array(Z)
        #artists[g, 1].set_array(N)

        Z = P - S
        #m = R < P
        #Z[m] = T[m] + S[m] - 2.0*P[m]
        artists[g, 2].set_array(Z)
        #artists[g, 2].set_array(N)

    if movie:
        fig.texts[2].set_text(f'low:{low:.2f}')
    return artists.flatten()

# Data

alphas = np.linspace(my.alphamax, my.alphamin, num=nr)
logess = np.linspace(my.logesmin, my.logesmax, num=nc)
rhos = 1.0 - 1.0/pow(2.0, logess)
RR, AA = np.meshgrid(rhos, alphas)
aBlows = np.linspace(0.1, 0.8, num=nframes)

# Figure properties

width = plotsize*len(titles)
height = plotsize*len(givens)
xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Influence of $\it{B}$'
biglabels = plotsize*5 + height/4
ticklabels = plotsize*4
xticks = [0, nc/2 - 0.5, nc - 1]
yticks = [0, nr/2 - 0.5, nr - 1]
xmin = logess[0]
xmax = logess[-1]
ymin = alphas[-1]
ymax = alphas[0]
xticklabels = [f'{xmin:2.0f}',
               f'{(xmin + xmax)/2.0:2.0f}',
               f'{xmax:2.0f}']
yticklabels = [f'{ymax:3.1f}',
               f'{(ymin + ymax)/2.0:3.1f}',
               f'{ymin:3.1f}']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Create figure

fig, axs = plt.subplots(nrows=len(givens),
                        ncols=len(titles),
                        figsize=(width, height))

left_x = axs[0, 0].get_position().x0
right_x = axs[-1, -1].get_position().x1
center_x = (left_x + right_x) / 2
top_y = axs[0, 0].get_position().y1
bottom_y = axs[-1, -1].get_position().y0
center_y = (top_y + bottom_y) / 2
fig.supxlabel(xlabel,
              x=center_x,
              y=bottom_y*0.3,
              fontsize=biglabels)
fig.supylabel(ylabel,
              x=left_x*0.4,
              y=center_y,
              fontsize=biglabels)

letterposition = 1.035
for i, ax in enumerate(fig.get_axes()):
    ax.set(xticks=xticks, yticks=yticks)
    ax.set(xticklabels=[], yticklabels=[])
    letter = ord('a') + i
    ax.text(0,
            letterposition,
            chr(letter),
            transform=ax.transAxes,
            fontsize=plotsize*5,
            weight='bold')
for g, given in enumerate(givens):
    axs[g, 0].set_yticklabels(yticklabels, fontsize=ticklabels)
for c, title in enumerate(titles):
    axs[0, c].set_title(title, pad=plotsize*10, fontsize=plotsize*5)
    axs[-1, c].set_xticklabels(xticklabels, fontsize=ticklabels)

if movie:
    fig.text(right_x,
             bottom_y*0.5,
             't\n0',
             fontsize=biglabels,
             color='grey',
             ha='right')

# Assign axs objects to variables
# (AxesImage)

artists = np.empty_like(axs) 
dummy_Z = np.empty_like(AA)
frames = aBlows
frame0 = aBlows[0]

for g, given in enumerate(givens):
    for c, title in enumerate(titles):
        artists[g, c] = axs[g, c].imshow(dummy_Z,
                                         vmin=0,
                                         vmax=1)

# Add data and save figure

if movie:
    ani = FuncAnimation(fig,
                        update,
                        frames=frames,
                        fargs=(artists,),
                        blit=True)
    ani.save(filename + '.mp4', writer='ffmpeg', fps=10)
else:
    update(frame0, artists,)
    plt.savefig(filename + '.png', transparent=False)

plt.close()

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
