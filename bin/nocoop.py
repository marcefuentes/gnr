#! /usr/bin/env python

from glob import glob
from math import log
import os
import sys
import time
import imageio.v2 as iio
import matplotlib.pyplot as plt
import numpy as np
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

letter = ('a', 'b', 'c', 'd', 'e', 'f')

dfglos = pd.DataFrame(
    [('ES', 'Substitutability of $\it{A}$', True, pow(2, -5.5), None), 
    ('alpha', 'Substitutability of $\it{A}$', True, None, None), 
    ('Given', 'Partner\'s share of $\it{A}$', False, None, None),
    ('ChooseCost', 'Cost of comparing potential partners', True, None, None),
    ('MimicCost', 'Cost of comparing partner to self', True, None, None),
    ('DeathRate', 'Death rate', True, 0.005, 0.2),
    ('GroupSize', 'Number of potential partners', True, 48.0, 2.6),
    ('MimicGrainInit', 'Sensitivity for mimicking partner', False, None, None),
    ('ChooseGrainInit', 'Sensitivity for choosing partner', False, None, None),
    ('N', 'Population size', True, None, None)],
    columns = ['name', 'label', 'log', 'min', 'max'])

dftraits = pd.DataFrame(
    [('ChooseGrainmedian', 'Sensitivity for\nchoosing partner', 600.0),
    ('ChooseGrain', 'Sensitivity for\nchoosing partner', 600.0),
    ('ChooseGrain12', 'Frequency of\npartner choosers', 600.0),
    ('MimicGrainmedian', 'Sensitivity for\nmimicking partner', 600.0),
    ('MimicGrain', 'Sensitivity for\nmimicking partner', 600.0),
    ('MimicGrain12', 'Frequency of\nreciprocators', 600.0),
    ('helpmedian', 'Help', 600.0*2.00),
    ('help', 'Help', 600.0*2.00),
    ('a2Seenmedian', 'Effort to get $\it{A}$', 600.0),
    ('a2Seen', 'Effort to get $\it{A}$', 600.0),
    ('a2Seen12', 'Frequency of\ndefectors', 600.0),
    ('a2Seen31', 'Frequency of\ncooperators', 600.0),
    ('a2Defaultmedian', 'Default $\it{a}$', 600.0),
    ('a2Default', 'Default $\it{a}$', 600.0),
    ('a2Default12', 'Frequency of\ndefectors', 600.0),
    ('a2Default31', 'Frequency of\ncooperators', 600.0),
    ('wmedian', 'Fitness', 600.0),
    ('w', 'Fitness', 600.0),
    ('chose_partner', 'Frequency of\nswitching to a new partner', 2000.0),
    ('changed_a2', 'Frequency of\nchanging $\it{a}$', 2000.0),
    ('helpBD', 'Fluctuation of help', 2000.0),
    ('wBD', 'Fluctuation of fitness', 2000.0)],
    columns = ['name', 'label', 'bubble_size'])

dfglos = dfglos.set_index('name')
dftraits = dftraits.set_index('name')

alpha = 0.5
R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0
b = a2max/a1max
npoints = 32

log_ess = np.linspace(-5, 5, num=11) if module.ic == 'ces' else np.linspace(-10, 0, num=11)
ess = pow(2, log_ess)
givens = np.linspace(1.0, 0.0, num=11)
givens[0] = 0.99999
x = np.linspace(1.0/(npoints + 1), 1.0 - 1.0/(npoints + 1), num = npoints)
a2partner = x
y = np.linspace(1.0 - 1.0/(npoints + 1), 1.0/(npoints + 1), num = npoints)
X, Y = np.meshgrid(x, y)

def fitness(x, y, given, rho):
    q1 = (1.0-y)*R1
    q2 = y*R2*(1.0-given) + x*R2*given
    if sys.argv[1] == 'q':
        w = 4.0*pow(q1, rho)/9.0 + 4.0*q2/9.0
    else:
        if rho == 0.0:
            w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
        else:
            w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

def a2maxw(given, rho):
    R = R2/R1
    T = b*R*(1.0 - given)
    Q = R*pow(T*(1.0 - alpha)/alpha, 1.0/(rho - 1.0))
    a1 = Q*(b - given + a2partner*given)/(1.0 + Q*b*(1 - given))
    a2 = a2max - b*a1
    a2 = 1.0 - a2
    return a2

class BarsAll:

    def prepare(self, dfs):

        bincount = int(sum(map(lambda x: module.top_traits[0] in x, [*dfs[folderlist[0]]]))/2) - 2

        self.traits = [{'name': trait} for trait in module.top_traits]
        for trait in self.traits:
            trait['namebins_list'] = [trait['name'] + str(x) for x in range(bincount)]
            trait['title'] = dftraits.loc[trait['name'], 'label']

        self.innercols = dfs[folderlist[0]][module.glos['x']].unique()
        self.innerrows = dfs[folderlist[0]][module.glos['y']].unique()
        self.innercols.sort()
        self.innerrows.sort() if module.glos['y'] == 'GroupSize' else self.innerrows[::-1].sort()

        self.bins = [(x+1)/bincount for x in range(bincount)]
        self.barwidth = 2.0/bincount

        self.color_blue = [red-0.95, green-0.95, blue-0.05]
        self.color_green = [red-0.95, green-0.05, blue-0.95]
        color_gray = [red-0.15, green-0.15, blue-0.15]
        colorsd_gray = [red-0.10, green-0.10, blue-0.10]

        self.colors = [color_gray, self.color_blue]
        self.alphas = [1.0, 0.9]

        return self

    def chart(self, dfts):
 
        fig = plt.figure(figsize=(width*2.0+1.2, height+0.5))
        hs=-0.5300
        ws=0.0
        fig.supxlabel(t=dfglos.loc[module.glos['x'], 'label'], x=0.51, y=0.00, fontsize=fslabel)
        fig.supylabel(t=dfglos.loc[module.glos['y'], 'label'], x=0.04, fontsize=fslabel, ha='center')

        if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=fstick, color='grey', ha='right')

        outergrids = fig.add_gridspec(nrows=1, ncols=2, wspace=0.1)

        innergrid = outergrids[0].subgridspec(nrows=len(self.innerrows), ncols=len(self.innercols), wspace=ws, hspace=hs)
        plt.subplots_adjust(hspace=0.0, wspace=0.0)
        axs = innergrid.subplots()
        axs[0, int(len(self.innercols)/2)].set_title('Fitness', pad=1.0, fontsize=fslabel) # Prints the title of the middle column. Bad if there are even columns
        
        for row, given in zip(axs, givens):
            for ax, es, log_es in zip(row, ess, log_ess):
                rho = 1.0 - 1.0/es if module.ic == 'ces' else es
                Z = fitness(X, Y, given, rho)
                ax.imshow(Z, vmin=0, vmax=2)
                xaxis = a2partner*npoints
                yaxis = a2maxw(given, rho)*npoints
                ax.plot(xaxis, yaxis, color='white')
                ax.set(xticks=[], yticks=[], xlim=(0, npoints-1), ylim=(npoints-1, 0))
                if given == 0.0:
                    ax.set_xlabel(round(log_es), fontsize=fstick)
                if log_es == -5:
                    ax.set_ylabel(round(given, 1), rotation='horizontal', horizontalalignment='right', verticalalignment='center', fontsize=fstick)
        fig.text(0.125, 0.86, 'a', fontsize=fslabel, weight='bold')

        innergrid = outergrids[1].subgridspec(nrows=len(self.innerrows), ncols=len(self.innercols), wspace=ws, hspace=hs)
        axs = innergrid.subplots()
        trait = self.traits[0]
        axs[0, int(len(self.innercols)/2)].set_title(trait['title'], pad=1.0, fontsize=fslabel) # Prints the title of the middle column. Bad if there are even columns
        ds = [dfts[module.folders[0][0]['control']], dfts[module.folders[0][0]['treatment']]]

        for row, (rowax, innerrow) in enumerate(zip(axs, self.innerrows)): 
            for column, (ax, innercol) in enumerate(zip(rowax, self.innercols)):
                medians = []
                [medians.append(d.loc[(d[module.glos['x']] == innercol) & (d[module.glos['y']] == innerrow), trait['name'] + 'median'].values[0]) for d in ds]
                dif = medians[0]/medians[1] + 0.05
                #dif = 0.5
                if 'Sensitivity' in trait['title']: dif = 1.0/dif
                self.colors[1] = dif_color(dif)
                for b, name0, name1 in zip(self.bins[::2], trait['namebins_list'][::2], trait['namebins_list'][1::2]):
                    barheight=ds[1].loc[(ds[1][module.glos['x']] == innercol) & (ds[1][module.glos['y']] == innerrow), name0] + ds[1].loc[(ds[1][module.glos['x']] == innercol) & (ds[1][module.glos['y']] == innerrow), name1]
                    ax.bar(x=b, height=barheight, align='edge', color=self.colors[1], linewidth=0, width=self.barwidth)
                    ax.set(xticks=[], yticks=[], xlim=[0, 1], ylim=[0, module.ylim])
                    ax.set_box_aspect(1)
                if row == len(self.innerrows) - 1:
                    xlabel = round(log(innercol, 2)) if dfglos.loc[module.glos['x'], 'log'] else innercol
                    ax.set_xlabel(xlabel, fontsize=fstick)
        fig.text(0.53, 0.86, 'b', fontsize=fslabel, weight='bold')

        plt.savefig(outfile, dpi=100)
        plt.close()

def dif_color(dif):
    color = (red*pow(dif,1.7), green*pow(dif,0.7), blue*pow(dif,1.7)) if dif <= 1.0 else (red*pow(dif,-1.7), green*pow(dif,-1.7), blue*pow(dif,-0.7))
    return color

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

if module.ftype == 'barsone':
    folderlist = folderlist_csv(folderlist)
    extension = '*.csv'
    pr = BarsOne()
elif module.ftype == 'barsall':
    folderlist = folderlist_csv(folderlist)
    extension = '*.csv'
    pr = BarsAll()
elif module.ftype == 'bubbles':
    folderlist = folderlist_csv(folderlist)
    extension = '*.csv'
    pr = Bubbles()
elif module.ftype == 'scatterall':
    folderlist = module.top_folders
    extension = '*.ics'
    pr = ScatterAll()
elif module.ftype == 'scatterone':
    folderlist = module.top_folders
    extension = '*.ics'
    pr = ScatterOne()
else:
    print('No such ftype')
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

if module.movie:
    frames = []
    for t in dfs[folderlist[0]].Time.unique():
        outfile = f'delete{t}.png'
        create_figure(t)
        frames.append(iio.imread(outfile))
        os.remove(outfile)
        percent = t*100/laststep
        print(f'Created {percent:.1f}% of frames', end='\r')
    print('\nAdding frames to movie...')
    iio.mimsave(f'{module.filename}.gif', frames)
else:
    outfile = f'{module.filename}.png'
    create_figure(laststep)

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
