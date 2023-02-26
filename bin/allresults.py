#! /usr/bin/env python

from glob import glob
import imageio.v2 as iio
import matplotlib.pyplot as plt
import mymodule
import numpy as np
import os
import pandas as pd
import time

start_time = time.perf_counter ()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

traits = ['a2Seenmean',
          'ChooseGrainmean',
          'MimicGrainmean',
          'wmean']
titles = ['Effort to get $\it{B}$',
               'Sensitivity for\nchoosing partner',
               'Sensitivity for\nmimicking partner',
               'Fitness']
traitvmaxs = [mymodule.a2max,
                mymodule.a2max,
                mymodule.a2max,
                mymodule.fitness(np.array([mymodule.a2max]),
                                    np.array([mymodule.a2max]),
                                    np.array([0.0]),
                                    np.array([0.9]),
                                    np.array([5.0]))]
folders = ['given0', 'none', 'p', 'r', 'pr', 'p8r']

movie = False
rows = folders
plotsize = 4

dfs = []
for folder in folders:
    df = pd.concat(map(pd.read_csv, glob(os.path.join(folder, '*.csv'))),
                    ignore_index=True)
    df.ChooseGrainmean = 1.0 - df.ChooseGrainmean
    df.MimicGrainmean = 1.0 - df.MimicGrainmean
    dfs.append(df)

df = dfs[1]
ts = df.Time.unique()
if movie:
    frames = []
else:
    ts = [ts[-1]]
rowindex = 'alpha'
nr = df['alpha'].nunique()
nc = df['logES'].nunique()

xlabel = 'Substitutability of $\it{B}$'
ylabel = 'Value of $\it{B}$'
letter = ord('a')
letterposition = 1.035
xmin = df['logES'].min()
xmax = df['logES'].max()
ymin = df['alpha'].min()
ymax = df['alpha'].max()
xticks = [-0.5, nc/2-0.5, nc-0.5]
yticks = [-0.5, nr/2-0.5, nr-0.5]
xticklabels = [f'{xmin:2.0f}',
               f'{(xmin + xmax)/2.0:2.0f}',
               f'{xmax:2.0f}']
yticklabels = [f'{ymax:3.1f}',
               f'{(ymin + ymax)/2.0:3.1f}',
               f'{ymin:3.1f}']
width = plotsize*len(titles)
height = plotsize*len(rows)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig, axs = plt.subplots(nrows=len(rows),
                        ncols=len(titles),
                        figsize=(width, height))
fig.supxlabel(xlabel,
              x=0.513,
              y=0.04,
              fontsize=width*3)
fig.supylabel(ylabel,
              x=0.03,
              y=0.493,
              fontsize=width*3)

for ax in fig.get_axes():
    ax.set(xticks=xticks,
           yticks=yticks,
           xticklabels=[],
           yticklabels=[])
    ax.text(0,
            letterposition,
            chr(letter),
            transform=ax.transAxes,
            fontsize=plotsize*5,
            weight='bold')
    letter += 1
for i, row in enumerate(rows):
    axs[i, 0].set_yticklabels(yticklabels, fontsize=plotsize*4)
for j, title in enumerate(titles):
    axs[0, j].set_title(title, pad=plotsize*10, fontsize=plotsize*6)
    axs[-1, j].set_xticklabels(xticklabels, fontsize=plotsize*4)

#        if letter <= ord('z'): 
#            textl = chr(letter)
#        else:
#            textl = 'a' + chr(letter - 26)
#        letter += 1

for t in ts:
    for i, df in enumerate(dfs):
        for j, trait in enumerate(traits):
            Z = pd.pivot_table(df.loc[df.Time == t],
                               values=trait,
                               index=[rowindex],
                               columns=['logES']).sort_index(axis=0,
                                                    ascending=False)
            axs[i, j].imshow(Z, vmin=0, vmax=traitvmaxs[j])
    if movie:
        text = fig.text(0.90,
                        0.90,
                        f't\n{t}',
                        fontsize=width*2,
                        color='grey',
                        ha='right')
        plt.savefig('temp.png', transparent=False)
        text.remove()
        frames.append(iio.imread('temp.png'))
        os.remove('temp.png')
    else:
        plt.savefig(filename + '.png', transparent=False)

plt.close()

if movie:
    iio.mimsave(filename + '.gif', frames)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
