#! /usr/bin/env python

from glob import glob
from math import log
import numpy as np
import os
import sys
import time
import imageio.v2 as iio
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
fslarge = 26 # Label font size
fssmall = 18 # Tick font size
red = 0.97
green = 0.97
blue = 0.97

letter = ('a', 'b', 'c', 'd', 'e', 'f')

dfglos = pd.DataFrame(
    [('ES', 'Substitutability of $\it{B}$', True, None), 
    ('alpha', 'Substitutability of $\it{B}$', True, None, None), 
    ('Given', 'Partner\'s share of $\it{B}$', False, None, None),
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
    ('a2Seenmedian', 'Effort to get $\it{B}$', 600.0),
    ('a2Seen', 'Effort to get $\it{B}$', 600.0),
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

class Bubbles:

    def prepare(self, dfs):

        self.glosx_min = None if pd.isnull(dfglos.loc[module.glos['x'], 'min']) else dfglos.loc[module.glos['x'], 'min']
        self.glosx_max = None if pd.isnull(dfglos.loc[module.glos['x'], 'max']) else dfglos.loc[module.glos['x'], 'max']
        self.glosy_min = None if pd.isnull(dfglos.loc[module.glos['y'], 'min']) else dfglos.loc[module.glos['y'], 'min']
        self.glosy_max = None if pd.isnull(dfglos.loc[module.glos['y'], 'max']) else dfglos.loc[module.glos['y'], 'max']

        return self

    def chart(self, dfts):

        fig, axs = plt.subplots(nrows=len(module.folders), ncols=len(module.top_traits), figsize=(width-1, height), sharex=True, sharey=True, constrained_layout=False, squeeze=False)
        fig.supxlabel(t=dfglos.loc[module.glos['x'], 'label'], y=0.02, x=0.513, fontsize=fslarge)
        fig.supylabel(t=dfglos.loc[module.glos['y'], 'label'], x=0.06, fontsize=fslarge, ha='center')

        if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=fssmall, color='grey', ha='right')

        [dfts[folder].sort_values(by=[module.glos['x'], module.glos['y']], inplace=True) for folder in folderlist]

        count = 0;
        for rowax, traits, folders in zip(axs, module.traits, module.folders):
            for ax, trait, folder in zip(rowax, traits, folders):
                df = dfts[folder['treatment']]
                x = df[module.glos['x']]
                if dfglos.loc[module.glos['x'], 'log']: x = x.apply(lambda i: log(i, 2))
                y = df[module.glos['y']]
                #if traits == module.top_traits: ax.set_title(dftraits.loc[trait, 'label'], pad=10.0, fontsize=fslarge)
                ax.set_title(dftraits.loc[trait, 'label'], pad=10.0, fontsize=fslarge)
                ax.tick_params(axis='x', labelsize=fssmall)
                ax.tick_params(axis='y', labelsize=fssmall)
                ax.set_xticks([-5, 0, 5])
                ax.set_xticklabels([-5, 0, 5])
                #if dfglos.loc[module.glos['x'], 'log']: ax.set_xscale('log', base=2)
                #if dfglos.loc[module.glos['y'], 'log']: ax.set_yscale('log', base=2)
                s = df[trait]
                difs = dfts[folder['control']][trait]/s
                if 'Sensitivity' in dftraits.loc[trait, 'label']: difs = 1.0/difs
                color = []
                [color.append(dif_color(dif)) for dif in difs]
                for suffix, suffixalpha in zip(['SD', ''], [0.2, 1.0]):
                    size = s + df[trait + suffix] if suffix == 'SD' else s
                    ax.scatter(x, y, c=color, ec=color, alpha=suffixalpha, s=size*dftraits.loc[trait, 'bubble_size'])
                    ax.set_xlim(self.glosx_min, self.glosx_max)
                    ax.set_ylim(self.glosy_min, self.glosy_max)
                    ax.set_box_aspect(1)
                ax.text(-6.6, 1.09, letter[count], fontsize=fslarge, weight='bold')
                count = count + 1

        plt.savefig(outfile, transparent=False)
        plt.close()

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
 
        fig = plt.figure(figsize=(width + 3.9, height + 1.0))
        fig.supxlabel(t=dfglos.loc[module.glos['x'], 'label'], y=0.00, fontsize=fslarge)
        fig.supylabel(t=dfglos.loc[module.glos['y'], 'label'], x=0.003*width, fontsize=fslarge, ha='center')

        if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=fssmall, color='grey', ha='right')

        outergrids = fig.add_gridspec(nrows=1, ncols=len(module.top_traits), wspace=0.1)

        ds = [dfts[module.folders[0][0]['control']], dfts[module.folders[0][0]['treatment']]]

        for trait, outergrid in zip(self.traits, outergrids):
            innergrid = outergrid.subgridspec(nrows=len(self.innerrows), ncols=len(self.innercols), wspace=0.0, hspace=0.0)
            axs = innergrid.subplots()
            axs[0, int(len(self.innercols)/2)].set_title(trait['title'], pad=1.0, fontsize=fslarge) # Prints the title of the middle column. Bad if there are even columns
            for row, (rowax, innerrow) in enumerate(zip(axs, self.innerrows)): 
                for column, (ax, innercol) in enumerate(zip(rowax, self.innercols)):
                    medians = []
                    [medians.append(d.loc[(d[module.glos['x']] == innercol) & (d[module.glos['y']] == innerrow), trait['name'] + 'median'].values[0]) for d in ds]
                    dif = medians[0]/medians[1]
                    if 'Sensitivity' in trait['title']: dif = 1.0/dif
                    self.colors[1] = dif_color(dif)
                    for d, color, alpha in zip(ds, self.colors, self.alphas):
                        for b, name0, name1 in zip(self.bins[::2], trait['namebins_list'][::2], trait['namebins_list'][1::2]):
                            barheight=d.loc[(d[module.glos['x']] == innercol) & (d[module.glos['y']] == innerrow), name0] + d.loc[(d[module.glos['x']] == innercol) & (d[module.glos['y']] == innerrow), name1]
                            ax.bar(x=b, height=barheight, align='edge', color=color, linewidth=0, width=self.barwidth, alpha=alpha)
                            ax.set(xticks=[], yticks=[], ylim=[0, module.ylim])
                            ax.set_box_aspect(1)
                    if (trait['name'] == module.top_traits[0]) & (column == 0):
                        y = '$2^{{{}}}$'.format(round(log(innerrow, 2))) if dfglos.loc[module.glos['y'], 'log'] else innerrow
                        ax.set_ylabel(y, rotation='horizontal', horizontalalignment='right', verticalalignment='center')
                    if row == len(self.innerrows) - 1:
                        x = '$2^{{{}}}$'.format(round(log(innercol, 2))) if dfglos.loc[module.glos['x'], 'log'] else innercol
                        ax.set_xlabel(x)

        plt.savefig(outfile, dpi=100)
        plt.close()

class BarsOne:

    def prepare(self, dfs):

        bincount = int(sum(map(lambda x: module.top_traits[0] in x, [*dfs[folderlist[0]]]))/2) - 2

        self.traits = [{'name': trait} for trait in module.top_traits]
        for trait in self.traits:
            trait['namebins_list'] = [trait['name'] + str(x) for x in range(bincount)]
            trait['namebinsds_list'] = [n + 'SD' for n in trait['namebins_list']]
            trait['max'] = 2.0 if trait['name'] == 'w' else 1.0 # For a1Max = a2Max = 1.0 & R1 = R2 = 2.0.
            trait['binslist'] = [(x+1)*trait['max']/bincount for x in range(bincount)]
            trait['barwidth'] = -trait['max']/bincount
            trait['label'] = dftraits.loc[trait['name'], 'label'] 

        return self

    def chart(self, dfts):

        fig, axs = plt.subplots(nrows=len(module.folders), ncols=len(module.top_traits), figsize=(width, height), sharey=True, constrained_layout=False, squeeze=False)
        fig.supylabel('Frequency', fontsize=fslarge, ha='center')
        if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=fssmall, color='grey', ha='right')

        for rowax, folders in zip(axs, module.folders):
            for ax, trait, folder in zip(rowax, self.traits, folders):
                if folders == module.bottom_folders: ax.set_xlabel(trait['label'], fontsize=fslarge)
                dif = dfts[folder['control']][trait['name'] + 'median'].values[0]/dfts[folder['treatment']][trait['name'] + 'median'].values[0]
                if 'Sensitivity' in trait['label']: dif = 1.0/dif
                color = dif_color(dif)
                for b, namebin, namebinsd in zip(trait['binslist'], trait['namebins_list'], trait['namebinsds_list']):
                    barheight = dfts[folder['treatment']][namebin]
                    barheightsd = dfts[folder['treatment']][namebinsd]
                    ax.bar(x=b, height=barheight, align='edge', color=color, linewidth=0.0, width=trait['barwidth'])
                    ax.bar(x=b, height=barheightsd, align='edge', color=color, linewidth=0.0, width=trait['barwidth'], bottom=barheight, alpha=0.2)
                ax.set(ylim=(0.0, module.ylim), yticks=(0.0, module.ylim), yticklabels=(0.0, module.ylim))
                ax.set(xlim=(0.0, trait['max']), xticks=(0.0, trait['max']), xticklabels=(0.0, trait['max']))
                ax.tick_params(axis='x', labelsize=fssmall)
                ax.tick_params(axis='y', labelsize=fssmall)
                if (len(module.folders) > 1) & (folders != module.bottom_folders): ax.set(xticks=[])
                ax.set_box_aspect(1)

        plt.savefig(outfile, dpi=100)
        plt.close()

class ScatterAll:

    def prepare(self, dfs):

        self.traits = module.top_traits
        for trait in self.traits:
            trait['xlimit'] = 2.0 if trait['x'] == 'w' else 1.0
            trait['ylimit'] = 2.0 if trait['y'] == 'w' else 1.0
            trait['title'] = '$\it{x}$ = ' + dftraits.loc[trait['x'], 'label'] + '\n$\it{y}$ = ' + dftraits.loc[trait['y'], 'label']

        self.innercols = dfs[folderlist[0]][module.glos['x']].unique()
        self.innerrows = dfs[folderlist[0]][module.glos['y']].unique()
        self.innercols.sort()
        self.innerrows.sort() if module.glos['y'] == 'GroupSize' else self.innerrows[::-1].sort()

        return self

    def chart(self, dfts):
 
        fig = plt.figure(figsize=(width + 0.0, height))
        fig.supxlabel(t=dfglos.loc[module.glos['x'], 'label'], y=0.00, fontsize=fslarge)
        fig.supylabel(t=dfglos.loc[module.glos['y'], 'label'], x=0.003*width, fontsize=fslarge, ha='center')
        if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=fssmall, color='grey', ha='right')
        outer_grid = fig.add_gridspec(nrows=len(module.folders), ncols=len(module.top_traits), hspace=0.1, wspace=0.1)

        for nr, row in enumerate(module.folders):
            dft = dfts[row]
            for nc, trait in enumerate(self.traits):
                innergrid = outer_grid[nr, nc].subgridspec(nrows=len(self.innerrows), ncols=len(self.innercols), wspace=0.0, hspace=0.0)
                axs = innergrid.subplots()
                if nr == 0: axs[0, int(len(self.innercols)/2)].set_title(trait['title'], fontsize=fslarge) # Prints the title of the middle column
                for row, (rowax, innerrow) in enumerate(zip(axs, self.innerrows)): 
                    for column, (ax, innercol) in enumerate(zip(rowax, self.innercols)):
                        x = dft.loc[(dft[module.glos['x']] == innercol) & (dft[module.glos['y']] == innerrow), trait['x']]
                        y = dft.loc[(dft[module.glos['x']] == innercol) & (dft[module.glos['y']] == innerrow), trait['y']]
                        alphas = dft.loc[(dft[module.glos['x']] == innercol) & (dft[module.glos['y']] == innerrow), module.zalpha]
                        alphas = 1.0-alphas if module.sensitive else alphas
                        ax.scatter(x, y, c='k', alpha=alphas, s=0.00001)
                        ax.set_xlim(0.0, trait['xlimit'])
                        ax.set_ylim(0.0, trait['ylimit'])
                        ax.tick_params(axis='x', labelsize=fssmall)
                        ax.tick_params(axis='y', labelsize=fssmall)
                        ax.set(xticks=[], yticks=[])
                        if (nc == 0) & (column == 0):
                            y = '$2^{{{}}}$'.format(round(log(innerrow, 2))) if dfglos.loc[module.glos['y'], 'log'] else innerrow
                            ax.set_ylabel(y, rotation='horizontal', horizontalalignment='right', verticalalignment='center')
                        if (nr == 1) & (row == len(self.innerrows) - 1):
                            x = '$2^{{{}}}$'.format(round(log(innercol, 2))) if dfglos.loc[module.glos['x'], 'log'] else innercol
                            ax.set_xlabel(x)
                        ax.set_box_aspect(1)

        plt.savefig(outfile, dpi=100)
        plt.close()

class ScatterOne:

    def prepare(self, dfs):

        self.traits = module.top_traits
        for trait in self.traits:
            trait['xlimit'] = 2.0 if trait['x'] == 'w' else 1.0
            trait['ylimit'] = 2.0 if trait['y'] == 'w' else 1.0
            trait['title'] = '$\it{x}$ = ' + dftraits.loc[trait['x'], 'label'] + '\n$\it{y}$ = ' + dftraits.loc[trait['y'], 'label']
        return self

    def chart(self, dfts):

        fig, axs = plt.subplots(nrows=len(module.folders), ncols=len(module.top_traits), figsize=(width, height), constrained_layout=False, squeeze=False)
        if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=fssmall, color='grey', ha='right')

        for rowax, row in zip(axs, module.folders):
            dft = dfts[row]
            for ax, trait in zip(rowax, self.traits):
                #ax.scatter(dft[trait['x']], dft[trait['y']], c=dft[module.zcolor], cmap=module.colormap, alpha=0.2, s=1)
                alphas = 1.0-dft[module.zalpha] if module.sensitive else dft[module.zalpha]
                ax.scatter(dft[trait['x']], dft[trait['y']], c='k', alpha=alphas, s=0.01)
                ax.set(xlim=(0.0, trait['xlimit']), xticks=(0.0, trait['xlimit']), xticklabels=(0.0, trait['xlimit']))
                ax.set(ylim=(0.0, trait['ylimit']), yticks=(0.0, trait['ylimit']), yticklabels=(0.0, trait['ylimit']))
                ax.tick_params(axis='x', labelsize=fssmall)
                ax.tick_params(axis='y', labelsize=fssmall)
                if row == module.top_row:    
                    ax.set_title(trait['title'], fontsize=fslarge)
                    ax.set(xticks=[])
                ax.set_box_aspect(1)

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