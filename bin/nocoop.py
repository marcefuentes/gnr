#! /usr/bin/env python

from glob import glob
from math import log
import numpy as np
import os
import sys
import time
import matplotlib.pyplot as plt
import pandas as pd

if len(sys.argv) < 2:
    print('You must add an argument with the python module name (for example, config)')
    exit(1)

start_time = time.perf_counter ()

sys.path.append(os.path.dirname(os.path.expanduser(sys.argv[1])))
module = __import__(sys.argv[1])

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

height = 6.0*len(module.folders)
width = 6.0*len(module.top_traits)
fslabel = 26 # Label font size
fstick = 18 # Tick font size
red = 0.97
green = 0.97
blue = 0.97

letters = ('a', 'b', 'c', 'd', 'e', 'f')

dfglos = pd.DataFrame(
    [('ES', 'Substitutability of $\it{A}$', True, None), 
    ('alpha', 'Substitutability of $\it{A}$', True, None, None), 
    ('Given', 'Partner\'s share of $\it{A}$', False, None, None),
    ('MimicGrainInit', 'Sensitivity for mimicking partner', False, None, None),
    ('ChooseGrainInit', 'Sensitivity for choosing partner', False, None, None),
    ('N', 'Population size', True, None, None)],
    columns = ['name', 'label', 'log', 'min', 'max'])

dftraits = pd.DataFrame(
    [('ChooseGrainmedian', 'Sensitivity for\nchoosing partner', 610.0, 1.0),
    ('ChooseGrain12', 'Frequency of\npartner choosers', 610.0),
    ('MimicGrainmedian', 'Sensitivity for\nmimicking partner', 610.0, 1.0),
    ('MimicGrain12', 'Frequency of\nreciprocators', 610.0),
    ('helpmedian', 'Help', 610.0, 1.0),
    ('a2Seenmedian', 'Effort to get $\it{A}$', 610.0, 0.5),
    ('a2Seen12', 'Frequency of\ndefectors', 610.0),
    ('a2Seen31', 'Frequency of\ncooperators', 610.0),
    ('wmedian', 'Fitness', 610.0, 1.0)],
    columns = ['name', 'label', 'bubble_size', 'vmax'])

dfglos = dfglos.set_index('name')
dftraits = dftraits.set_index('name')

alpha = 0.5
R1 = 2.0
R2 = 2.0
R = R2/R1
a1max = 1.0
a2max = 1.0
b = a2max/a1max

num = 2000
log_ess = np.linspace(-5, 5, num=num)
rhos = 1.0 - 1.0/pow(2, log_ess)
givens = np.linspace(1.0, 0.0, num=num)
givens[0] = 0.99999
Xrhos, Ygivens = np.meshgrid(rhos, givens)

class Bubbles:

    def prepare(self, dfs):

        self.glosx_min = None if pd.isnull(dfglos.loc[module.glos['x'], 'min']) else dfglos.loc[module.glos['x'], 'min']
        self.glosx_max = None if pd.isnull(dfglos.loc[module.glos['x'], 'max']) else dfglos.loc[module.glos['x'], 'max']
        self.glosy_min = None if pd.isnull(dfglos.loc[module.glos['y'], 'min']) else dfglos.loc[module.glos['y'], 'min']
        self.glosy_max = None if pd.isnull(dfglos.loc[module.glos['y'], 'max']) else dfglos.loc[module.glos['y'], 'max']

        return self

    def chart(self, dfts):

        fig, axs = plt.subplots(nrows=2, ncols=len(module.top_traits), figsize=(width-1, height*2), constrained_layout=False, squeeze=False)
        fig.supxlabel(t=dfglos.loc[module.glos['x'], 'label'], y=0.02, x=0.513, fontsize=fslabel)
        fig.supylabel(t=dfglos.loc[module.glos['y'], 'label'], x=0.06, fontsize=fslabel, ha='center')

        [dfts[folder].sort_values(by=[module.glos['x'], module.glos['y']], inplace=True) for folder in folderlist]

        extent = 0, num, 0, num
        Z0 = a2eq(Xrhos, Ygivens)
        axs[0, 0].imshow(Z0, extent=extent, cmap='magma', vmin=0, vmax=0.5)

        Z1 = Z0*R2*Ygivens
        axs[0, 1].imshow(Z1, extent=extent, cmap='magma', vmin=0, vmax=1.0)

        Z1 = fitness(Z0, Xrhos, Ygivens)
        axs[0, 2].imshow(Z1, extent=extent, cmap='magma', vmin=0, vmax=1.0)

        for ax, trait, folder in zip(axs[1], module.top_traits, module.top_folders):
            df = dfts[folder['treatment']]
            x = df[module.glos['x']]
            if dfglos.loc[module.glos['x'], 'log']: x = x.apply(lambda i: log(i, 2))
            y = df[module.glos['y']]
            s = df[trait]
            if 'Sensitivity' in dftraits.loc[trait, 'label']: s = 1 - s
            ax.scatter(x, y, c=s, cmap='magma', marker='s', vmin=0, vmax=dftraits.loc[trait, 'vmax'], s=dftraits.loc[trait, 'bubble_size'])
            ax.set_xlim(self.glosx_min, self.glosx_max)
            ax.set_ylim(self.glosy_min, self.glosy_max)
            ax.set_box_aspect(1)

        for ax, trait in zip(axs[0], module.top_traits):
            ax.set_title(dftraits.loc[trait, 'label'], pad=10.0, fontsize=fslabel)
            ax.set_xticks([0, num/2.0, num])
            ax.set_yticks([0, num/2.0, num])
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        for ax in axs[1]:
            ax.set_xticks([-5, 0, 5])
            ax.set_yticks([0.0, 0.5, 1.0])
            ax.set_xticklabels([-5, 0, 5], fontsize=fstick)
            ax.set_yticklabels([])
        axs[0, 0].set_yticklabels([0.0, 0.5, 1.0], fontsize=fstick) 
        axs[1, 0].set_yticklabels([0.0, 0.5, 1.0], fontsize=fstick) 
        for ax, letter in zip(axs[0], letters[:3]):
            ax.text(-200, 2070, letter, fontsize=fslabel, weight='bold')
        for ax, letter in zip(axs[1], letters[3:]):
            ax.text(-6.6, 1.09, letter, fontsize=fslabel, weight='bold')
            
        plt.savefig(outfile, transparent=False)
        plt.close()

def a2eq(X, Y):
    T = b*R*(1.0 - Y)
    Q = R*pow(T*(1.0 - alpha)/alpha, 1.0/(X - 1.0))
    a2 = 1.0/(1.0 + Q)
    return a2

def fitness(Z, X, Y):
    q1 = (1.0 - Z)*R1
    q2 = Z*R2*(1.0 - Y) + Z*R2*Y
    w = np.where(X == 0.0, pow(q1, alpha)*pow(q2, 1.0 - alpha), pow(alpha*pow(q1, X) + (1.0 - alpha)*pow(q2, X), 1.0/X)) 
    return w

def create_figure(t):
    dfts = {}
    for folder in folderlist:
        dfts[folder] = dfs[folder].loc[dfs[folder]['Time'] == t].copy()
    pr.chart(dfts)

def folderlist_csv(folderlist):
    for folders in module.folders:
        for folder in folders:
            folderlist.append(folder['treatment'])
            folderlist.append(folder['control'])
    folderlist = list(set(folderlist))
    return folderlist

folderlist = []

if module.ftype == 'bubbles':
    folderlist = folderlist_csv(folderlist)
    extension = '*.csv'
    pr = Bubbles()
else:
    print('ftype must be bubbles')
    exit(1)

dfs = {}
for folder in folderlist:
    dfs[folder] = pd.concat(map(pd.read_csv, glob(os.path.join(folder, extension))), ignore_index=True)

if 'one' in module.ftype:
    glovalue_x = float(str(f"{pow(2, int(module.glovalue['x'])):.6f}")) if dfglos.loc[module.glos['x'], 'log'] else module.glovalue['x'] 
    glovalue_y = float(str(f"{pow(2, int(module.glovalue['y'])):.6f}")) if dfglos.loc[module.glos['y'], 'log'] else module.glovalue['y']
    for folder in folderlist:
        dfs[folder] = dfs[folder].loc[(dfs[folder][module.glos['x']] == glovalue_x) & (dfs[folder][module.glos['y']] == glovalue_y)]

pr.prepare(dfs)

laststep = dfs[folderlist[0]].Time.iat[-1]

outfile = f'{module.filename}.png'
create_figure(laststep)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
