#!/usr/bin/python

import glob
import imageio
import math
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
module = __import__(sys.argv[1].replace('.py', ''))

if len(sys.argv) < 4:
    print('You must enter three arguments: scattergrain/barsgrain/fluct/change/fgrain/ffluct/fchange snap/movie uni/multi')
    exit(1)

fs = 16 # Label font size
lc = 0.9 # Lightest color
red = 1.0
green = 1.0
blue = 1.0

if sys.argv[3] == 'uni':
    outfile = f'{sys.argv[1]}_{module.choosecost}_{module.mimiccost}.png'
    choosecost = float(str("{:.6f}".format(pow(2, -int(module.choosecost)))))
    mimiccost = float(str("{:.6f}".format(pow(2, -int(module.mimiccost)))))
else:
    outfile = f'{sys.argv[1]}.png'

def create_figure(dfs, t):
    fig = plt.figure(figsize=(module.width, module.height))
    if sys.argv[3] == 'multi':
        fig.supxlabel(module.x_label, fontsize=fs)
        fig.supylabel(t=module.y_label, x=0.003*module.width, fontsize=fs)
    else:
        fig.supylabel(t='Frequency', x=0.003*module.width, fontsize=fs)
    if sys.argv[2] == 'movie':
        fig.text(0.89, 0.95, f'Time = {t}', fontsize=fs, ha='right')
    return fig

class barpr:

    inner_rows = []
    inner_cols = []
    xtick_label_lists = []
    c_name_suffixes = []
    c_names = []
    c_names_sd = []

    def prepare(self, dfs):

        self.bincount = int(sum(map(lambda x: module.c_name_roots[0] in x, [*dfs[0]]))/2) - 2
        self.c_name_suffixes = [x for x in range(self.bincount)]
        for root in module.c_name_roots:
            rootlist = []
            rootlist_sd = []
            for suffix in self.c_name_suffixes:
                rootlist.append(root + str(suffix)) 
                rootlist_sd.append(root + 'SD' + str(suffix)) 
            self.c_names.append(rootlist)
            self.c_names_sd.append(rootlist_sd)

        self.inner_cols = dfs[0][module.x_axis].unique()
        self.inner_rows = dfs[0][module.y_axis].unique()
        #self.inner_cols[::-1].sort()
        self.inner_cols.sort()
        self.inner_rows[::-1].sort()

        self.bins = [(x+1)/self.bincount for x in range(self.bincount)]
        self.bins_f = [(x+1)*2.0/self.bincount for x in range(self.bincount)]

        self.xticks = [0, 1]
        self.xtick_labels = [0, 1]
        self.xtick_labels_f = [0, 1.87]

        for color in module.colors:
            if color == 'grain':
                self.xtick_label_lists.append(self.xtick_labels)
            elif color == 'help':
                self.xtick_label_lists.append(self.xtick_labels)
            else:
                self.xtick_label_lists.append(self.xtick_labels_f)

        return self

    def chart(self, dfs, t):
 
        fig = create_figure(dfs, t)

        outer_grid = fig.add_gridspec(nrows=1, ncols=len(module.c_name_roots), wspace=0.1)

        width = 2.0/self.bincount

        for n, (c_name, ymax, title) in enumerate(zip(self.c_names, module.ymax, module.titles)):
            inner_grid = outer_grid[n].subgridspec(nrows=len(self.inner_rows), ncols=len(self.inner_cols), wspace=0, hspace=0)
            axs = inner_grid.subplots()
            #axs[int(len(self.inner_cols)/2)].set_title(title, fontsize=fs) # Prints the title of the middle column. Bad if there are even columns
            for row, (rowax, inner_row) in enumerate(zip(axs, self.inner_rows)): 
                for column, (ax, inner_col) in enumerate(zip(rowax, self.inner_cols)):
                    for b, name0, name1 in zip(self.bins[::2], c_name[::2], c_name[1::2]):
                        df = dfs[0]
                        height=df.loc[(df[module.x_axis] == inner_col) & (df[module.y_axis] == inner_row) & (df.Time == t), name0] + df.loc[(df[module.x_axis] == inner_col) & (df[module.y_axis] == inner_row) & (df.Time == t), name1]
                        ax.bar(x=b, height=height, align = 'edge', color=(red-0.15, green-0.15, blue-0.15), linewidth=0, width=width)
                        ax.set(xticks=[], yticks=[], ylim=[0, ymax])
                    for b, name0, name1 in zip(self.bins[::2], c_name[::2], c_name[1::2]):
                        df = dfs[1]
                        height=df.loc[(df[module.x_axis] == inner_col) & (df[module.y_axis] == inner_row) & (df.Time == t), name0] + df.loc[(df[module.x_axis] == inner_col) & (df[module.y_axis] == inner_row) & (df.Time == t), name1]
                        ax.bar(x=b, height=height, align = 'edge', color=(red-0.2, green-0.4, blue-0.7), linewidth=0, width=width)
                        ax.set(xticks=[], yticks=[], ylim=[0, ymax])
                    if (n == 0) & (column == 0):
                        if module.log == True:
                            y = '$2^{{{}}}$'.format(round(math.log(inner_row, 2)))
                        else:
                            y = inner_row
                        ax.set_ylabel(y, rotation='horizontal', horizontalalignment='right', verticalalignment='center')
                    if row == len(self.inner_cols) - 1:
                        if module.log == True:
                            x = '$2^{{{}}}$'.format(round(math.log(inner_col, 2)))
                        else:
                            x = inner_col
                        ax.set_xlabel(x)
                    if (column == len(self.inner_cols)/2):
                            ax.set_title(title, fontsize=fs)

        plt.savefig(outfile, dpi=100)
        plt.close()

    def uni(self, dfs, t, choosecost, mimiccost):

        fig = create_figure(dfs, t)

        width = -1.0/self.bincount

        outer_grid = fig.add_gridspec(nrows=1, ncols=len(module.c_name_roots), wspace=0.1)

        axs = outer_grid.subplots()

        for ax, title in zip(axs, module.titles):
            ax.set_xlabel(title, fontsize=fs)

        for column, (ax, name_root, c_name, c_name_sd, ymax, xtick_label_list) in enumerate(zip(axs, module.c_name_roots, self.c_names, self.c_names_sd, module.ymax, self.xtick_label_lists)):
            df = dfs[0]
            median0 = df.loc[(df[module.x_axis] == mimiccost) & (df[module.y_axis] == choosecost) & (df.Time == t), name_root + 'median'].values[0]
            df = dfs[1]
            median1 = df.loc[(df[module.x_axis] == mimiccost) & (df[module.y_axis] == choosecost) & (df.Time == t), name_root + 'median'].values[0]
            dif = median0 - median1
            if (name_root == 'ChooseGrain') or (name_root == 'MimicGrain') or ('BD' in name_root):
                dif = -dif
            if dif < 0.0:
                color=(red-0.95, green-0.95, blue-0.05)
                colorsd=(red-0.30, green-0.30, blue-0.05)
            else:
                color=(red-0.95, green-0.05, blue-0.95)
                colorsd=(red-0.30, green-0.05, blue-0.30)
            for b, name, namesd in zip(self.bins, c_name, c_name_sd):
                df = dfs[0]
                height = df.loc[(df[module.x_axis] == mimiccost) & (df[module.y_axis] == choosecost) & (df.Time == t), name]
                heightsd = df.loc[(df[module.x_axis] == mimiccost) & (df[module.y_axis] == choosecost) & (df.Time == t), namesd]
                ax.bar(x=b, height=height, align='edge', color=color, linewidth=0, width=width, alpha=1.0)
                ax.bar(x=b, height=heightsd, align='edge', color=colorsd, linewidth=0, width=width, bottom=height, alpha=1.0)
            ax.set(ylim=[0, ymax])
            ax.set(xticks=self.xticks, xticklabels=xtick_label_list)
            if (column > 0):
                ax.set(yticks=[])
            for b, name, namesd in zip(self.bins, c_name, c_name_sd):
                df = dfs[1]
                height = df.loc[(df[module.x_axis] == mimiccost) & (df[module.y_axis] == choosecost) & (df.Time == t), name]
                heightsd = df.loc[(df[module.x_axis] == mimiccost) & (df[module.y_axis] == choosecost) & (df.Time == t), namesd]
                ax.bar(x=b, height=height, align='edge', color=(red-0.15, green-0.15, blue-0.15), linewidth=0, width=width, alpha=0.9)
                ax.bar(x=b, height=heightsd, align='edge', color=(red-0.10, green-0.10, blue-0.10), linewidth=0, width=width, bottom=height, alpha=0.9)
            ax.set(ylim=[0, ymax])
            ax.set(xticks=self.xticks, xticklabels=xtick_label_list)
            if (column > 0):
                ax.set(yticks=[])

        plt.savefig(outfile, dpi=100)
        plt.close()

class scatterpr:

    def prepare(self, dfs):
        pass

    def chart(self, dfs, t):

        fig = create_figure(dfs, t)

        outer_grid = fig.add_gridspec(nrows=len(dfs), ncols=len(module.c_name_roots), wspace=0.1, hspace=0.1)

        axs = outer_grid.subplots(sharex='row', sharey='col')

        for ax, title in zip(axs[0], module.titles):
            ax.set_title(title, fontsize=fs)

        for row, (rowax, df) in enumerate(zip(axs, dfs)):
            for column, (ax, name_root) in enumerate(zip(rowax, module.c_name_roots)):
                if module.log == True:
                    ax.set_xscale('log', base=2)
                    ax.set_yscale('log', base=2)
                ax.set_aspect('equal', 'box')
                dif = (df.loc[df.Time == t, name_root] - dfs[1].loc[df.Time == t, name_root])
                if (name_root == 'ChooseGrainmedian') or (name_root == 'MimicGrainmedian') or ('BD' in name_root):
                    dif = -dif
                color = []
                for i in dif:
                    if i < 0.0:
                        color.append((red + i, green + i, blue))
                    else:
                        color.append((red - i, green, blue - i))
                x = df.loc[df.Time == t, module.x_axis]
                y = df.loc[df.Time == t, module.y_axis]
                s = df.loc[df.Time == t, name_root]
                ax.scatter(x, y, c=color, edgecolor='0.200', alpha=1.0, s=s*module.bubble_size)
                if (row == 0):
                    ax.set(xticks=[])
                if (column > 0):
                    ax.set(yticks=[])
                ax.set_xlim(module.x_min, module.x_max)
                ax.set_ylim(module.y_min, module.y_max)
                s = s + df.loc[df.Time == t, name_root + 'SD']
                ax.scatter(x, y, c=color, edgecolor='0.400', alpha=0.2, s=s*module.bubble_size)
                if (row == 0):
                    ax.set(xticks=[])
                if (column > 0):
                    ax.set(yticks=[])
                ax.set_xlim(module.x_min, module.x_max)
                ax.set_ylim(module.y_min, module.y_max)

        plt.savefig(outfile)
        plt.close()

def get_data(dfs):

    dirs = ['', 'drift/']

    for d in dirs:
        df = pd.concat(map(pd.read_csv, glob.glob(os.path.join(d, '*.csv'))), ignore_index=True)
        if sys.argv[2] != 'movie':
            lastt =  df.Time.iat[-1]
            df = df[df.Time == lastt]                
        dfs.append(df)

    return dfs

dfs = []
dfs = get_data(dfs)

if module.ftype == 'bars':
    pr = barpr()
else:
    pr = scatterpr()

pr.prepare(dfs)

if sys.argv[2] == 'movie':
    outfiles = []
    for t in dfs[0].Time.unique():
        print(f'Processing time {t}', end='\r')
        outfile = f'delete{t}.png'
        if sys.argv[3] == 'uni':
            pr.uni(dfs, t, choosecost, mimiccost)
        else:
            pr.chart(dfs, t)
        outfiles.append(outfile)
    if sys.argv[3] == 'multi':
        giffile = f'{sys.argv[1]}.gif'
    else:
        giffile = f'{sys.argv[1]}{module.choosecost}{module.mimiccost}.gif'
    with imageio.get_writer(giffile, mode='I') as writer:
        for outfile in outfiles:
            print(f'Adding {outfile} to movie', end='\r')
            image = imageio.imread(outfile)
            writer.append_data(image)
    for outfile in set(outfiles):
        os.remove(outfile)
else:
    if sys.argv[3] == 'multi':
        pr.chart(dfs, dfs[0].Time.iat[-1])
    else:
        pr.uni(dfs, dfs[0].Time.iat[-1], choosecost, mimiccost)

