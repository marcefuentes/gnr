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

if len(sys.argv) < 2:
    print('You must add an argument with the python module (for example, scattergrain)')
    exit(1)

outfile = f'{sys.argv[1]}.png'

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fslabel = 18 # Label font size
fstitle= 18 # Label font size
fstick = 14 # Label font size
red = 0.97
green = 0.97
blue = 0.97

class barsallpr:

    def prepare(self, dfs):

        self.bincount = int(sum(map(lambda x: module.c_name_roots[0] in x, [*dfs[0]]))/2) - 2

        self.c_name_suffixes = [x for x in range(self.bincount)]
        self.c_names = []
        for root in module.c_name_roots:
            rootlist = []
            for suffix in self.c_name_suffixes:
                rootlist.append(root + str(suffix)) 
            self.c_names.append(rootlist)

        self.inner_cols = dfs[0][module.x_axis].unique()
        self.inner_rows = dfs[0][module.y_axis].unique()
        self.inner_cols.sort()
        if module.y_axis == 'GroupSize':
            self.inner_rows.sort()
        else:
            self.inner_rows[::-1].sort()

        self.bins = [(x+1)/self.bincount for x in range(self.bincount)]

        return self

    def chart(self, dfs, t):
 
        fig = plt.figure(figsize=(module.width, module.height))
        fig.supxlabel(module.x_label, fontsize=fslabel)
        fig.supylabel(t=module.y_label, x=0.003*module.width, fontsize=fslabel)

        if module.movie == True:
            fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        outer_grid = fig.add_gridspec(nrows=1, ncols=len(module.c_name_roots), wspace=0.1)

        width = 2.0/self.bincount

        for n, (c_name, name_root, title) in enumerate(zip(self.c_names, module.c_name_roots, module.titles)):
            inner_grid = outer_grid[n].subgridspec(nrows=len(self.inner_rows), ncols=len(self.inner_cols), wspace=0, hspace=0)
            axs = inner_grid.subplots()
            axs[0, int(len(self.inner_cols)/2)].set_title(title, fontsize=fstick) # Prints the title of the middle column. Bad if there are even columns
            for row, (rowax, inner_row) in enumerate(zip(axs, self.inner_rows)): 
                for column, (ax, inner_col) in enumerate(zip(rowax, self.inner_cols)):
                    df = dfs[0]
                    median0 = df.loc[(df[module.x_axis] == inner_col) & (df[module.y_axis] == inner_row) & (df.Time == t), name_root + 'median'].values[0]
                    df = dfs[1]
                    median1 = df.loc[(df[module.x_axis] == inner_col) & (df[module.y_axis] == inner_row) & (df.Time == t), name_root + 'median'].values[0]
                    dif = median0 - median1
                    if (name_root == 'ChooseGrain') or (name_root == 'MimicGrain') or ('BD' in name_root):
                        dif = -dif
                    if dif < 0.0:
                        color=(red-0.95, green-0.95, blue-0.05)
                    else:
                        color=(red-0.95, green-0.05, blue-0.95)
                    for b, name0, name1 in zip(self.bins[::2], c_name[::2], c_name[1::2]):
                        df = dfs[0]
                        height=df.loc[(df[module.x_axis] == inner_col) & (df[module.y_axis] == inner_row) & (df.Time == t), name0] + df.loc[(df[module.x_axis] == inner_col) & (df[module.y_axis] == inner_row) & (df.Time == t), name1]
                        ax.bar(x=b, height=height, align='edge', color=color, linewidth=0, width=width)
                        ax.set(xticks=[], yticks=[], ylim=[0, module.ymax*2])
                    for b, name0, name1 in zip(self.bins[::2], c_name[::2], c_name[1::2]):
                        df = dfs[1]
                        height=df.loc[(df[module.x_axis] == inner_col) & (df[module.y_axis] == inner_row) & (df.Time == t), name0] + df.loc[(df[module.x_axis] == inner_col) & (df[module.y_axis] == inner_row) & (df.Time == t), name1]
                        ax.bar(x=b, height=height, align='edge', color=(red-0.15, green-0.15, blue-0.15), linewidth=0, width=width, alpha=0.9)
                        ax.set(xticks=[], yticks=[], ylim=[0, module.ymax*2])
                    if (n == 0) & (column == 0):
                        if module.logy == True:
                            y = '$2^{{{}}}$'.format(round(math.log(inner_row, 2)))
                        else:
                            y = inner_row
                        ax.set_ylabel(y, rotation='horizontal', horizontalalignment='right', verticalalignment='center')
                    if row == len(self.inner_rows) - 1:
                        if module.logx == True:
                            x = '$2^{{{}}}$'.format(round(math.log(inner_col, 2)))
                        else:
                            x = inner_col
                        ax.set_xlabel(x)
                    if column == len(self.inner_cols)/2:
                        ax.set_title(title, fontsize=fstitle)

        plt.savefig(outfile, dpi=100)
        plt.close()

class barsonepr:

    def prepare(self, dfs):

        self.x_value = float(str("{:.6f}".format(pow(2, int(module.x_value)))))
        self.y_value = float(str("{:.6f}".format(pow(2, int(module.y_value)))))

        self.bincount = int(sum(map(lambda x: module.c_name_roots[0] in x, [*dfs[0]]))/2) - 2

        self.c_name_suffixes = [x for x in range(self.bincount)]
        self.c_names = []
        self.c_names_sd = []
        for root in module.c_name_roots:
            rootlist = []
            rootlist_sd = []
            for suffix in self.c_name_suffixes:
                rootlist.append(root + str(suffix)) 
                rootlist_sd.append(root + 'SD' + str(suffix)) 
            self.c_names.append(rootlist)
            self.c_names_sd.append(rootlist_sd)

        self.bins = [(x+1)/self.bincount for x in range(self.bincount)]

        self.xticks = (0, 1)
        self.xtick_labels = (0, 1)
        self.xtick_labels_f = (0, 1.87)

        self.xtick_label_lists = []

        for color in module.colors:
            if color == 'grain':
                self.xtick_label_lists.append(self.xtick_labels)
            elif color == 'help':
                self.xtick_label_lists.append(self.xtick_labels)
            else:
                self.xtick_label_lists.append(self.xtick_labels_f)

        return self

    def chart(self, dfs, t):

        fig, axs = plt.subplots(nrows=1, ncols=len(module.c_name_roots), figsize=(module.width, module.height), constrained_layout=True)

        fig.supylabel('Frequency', fontsize=fslabel)

        if module.movie == True:
            fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        width = -1.0/self.bincount

        for ax, title in zip(axs, module.titles):
            ax.set_xlabel(title, fontsize=fslabel)

        for ax, name_root, c_name, c_name_sd, xtick_label_list in zip(axs, module.c_name_roots, self.c_names, self.c_names_sd, self.xtick_label_lists):
            df = dfs[0]
            median0 = df.loc[(df[module.x_axis] == self.x_value) & (df[module.y_axis] == self.y_value) & (df.Time == t), name_root + 'median'].values[0]
            df = dfs[1]
            median1 = df.loc[(df[module.x_axis] == self.x_value) & (df[module.y_axis] == self.y_value) & (df.Time == t), name_root + 'median'].values[0]
            dif = median0 - median1
            if (name_root == 'ChooseGrain') or (name_root == 'MimicGrain') or ('BD' in name_root):
                dif = -dif
            if dif < 0.0:
                color=(red-0.95, green-0.95, blue-0.05)
                colorsd=(red-0.30, green-0.30, blue-0.05)
            else:
                color=(red-0.95, green-0.05, blue-0.95)
                colorsd=(red-0.30, green-0.05, blue-0.30)
            alpha=1.0
            df = dfs[0]
            for b, name, namesd in zip(self.bins, c_name, c_name_sd):
                height = df.loc[(df[module.x_axis] == self.x_value) & (df[module.y_axis] == self.y_value) & (df.Time == t), name]
                heightsd = df.loc[(df[module.x_axis] == self.x_value) & (df[module.y_axis] == self.y_value) & (df.Time == t), namesd]
                ax.bar(x=b, height=height, align='edge', color=color, linewidth=0, width=width, alpha=alpha)
                ax.bar(x=b, height=heightsd, align='edge', color=colorsd, linewidth=0, width=width, bottom=height, alpha=alpha)
            ax.set(ylim=(0, module.ymax), yticks=(0, module.ymax), yticklabels=(0, module.ymax))
            ax.set(xticks=self.xticks, xticklabels=xtick_label_list)
            if name_root != module.c_name_roots[0]:
                ax.set(yticks=[])
            color=(red-0.15, green-0.15, blue-0.15)
            colorsd=(red-0.10, green-0.10, blue-0.10)
            alpha=0.9
            df = dfs[1]
            for b, name, namesd in zip(self.bins, c_name, c_name_sd):
                height = df.loc[(df[module.x_axis] == self.x_value) & (df[module.y_axis] == self.y_value) & (df.Time == t), name]
                heightsd = df.loc[(df[module.x_axis] == self.x_value) & (df[module.y_axis] == self.y_value) & (df.Time == t), namesd]
                ax.bar(x=b, height=height, align='edge', color=color, linewidth=0, width=width, alpha=alpha)
                ax.bar(x=b, height=heightsd, align='edge', color=colorsd, linewidth=0, width=width, bottom=height, alpha=alpha)
            ax.set(ylim=(0, module.ymax), yticks=(0, module.ymax), yticklabels=(0, module.ymax))
            ax.set(xticks=self.xticks, xticklabels=xtick_label_list)
            if name_root != module.c_name_roots[0]:
                ax.set(yticks=[])

        plt.savefig(outfile, dpi=100)
        plt.close()

class scatterpr:

    def prepare(self, dfs):
        self.suffixes = ('SD', '')
        self.suffixalphas = (0.2, 1.0)
        self.suffixecs = ('0.300', '0.000')

    def chart(self, dfs, t):

        fig, axs = plt.subplots(nrows=1, ncols=len(module.c_name_roots), figsize=(module.width, module.height), constrained_layout=False)

        fig.supxlabel(module.x_label, fontsize=fslabel)
        fig.supylabel(module.y_label, fontsize=fslabel, ha='center')

        if module.movie == True:
            fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        for ax, name_root, title in zip(axs, module.c_name_roots, module.titles):
            ax.set_title(title, fontsize=fstitle)
            ax.tick_params(axis='x', labelsize=fstick)
            ax.tick_params(axis='y', labelsize=fstick)
            if module.logx == True:
                ax.set_xscale('log', base=2)
            if module.logy == True:
                ax.set_yscale('log', base=2)
            dif = 1.0 - dfs[0].loc[dfs[0].Time == t, name_root]/dfs[1].loc[dfs[1].Time == t, name_root]
            if (name_root == 'ChooseGrainmedian') or (name_root == 'MimicGrainmedian') or ('BD' in name_root):
                dif = -dif
            color = []
            for i in dif:
                if i < 0.0:
                    if i < -red:
                        i = -red
                    color.append((red + i, green + i/2.0, blue + i))
                else:
                    if i > red:
                        i = red
                    color.append((red - i, green - i, blue - i/2.0))
            for drift, df in enumerate(dfs):
                for suffix, suffixec, suffixalpha in zip(self.suffixes, self.suffixecs, self.suffixalphas):
                    x = df.loc[df.Time == t, module.x_axis]
                    y = df.loc[df.Time == t, module.y_axis]
                    s = df.loc[df.Time == t, name_root]
                    if suffix == 'SD':
                        s = s + df.loc[df.Time == t, name_root + suffix]
                    if drift == 0:
                        ax.scatter(x, y, c=color, ec=color, alpha=suffixalpha, s=s*module.bubble_size)
                    else:
                        ax.scatter(x, y, facecolors='none', edgecolor=suffixec, alpha=0.05, s=s*module.bubble_size)
                    if name_root != module.c_name_roots[0]:
                        ax.set(yticks=[])
                    ax.set_xlim(module.x_min, module.x_max)
                    ax.set_ylim(module.y_min, module.y_max)
                    ax.set_box_aspect(1)

        plt.savefig(outfile, transparent=False)
        plt.close()

    def chart2(self, dfs, t):

        fig, axs = plt.subplots(nrows=2, ncols=len(module.c_name_roots), figsize=(module.width, module.height), constrained_layout=True)

        fig.supxlabel(module.x_label, fontsize=fslabel)
        fig.supylabel(module.y_label, fontsize=fslabel)

        if module.movie == True:
            fig.text(0.50, 0.50, f'Time = {t}', fontsize=fstitle, ha='center')

        for ax, title in zip(axs[0], module.titles):
            ax.set_title(title, fontsize=fstitle)

        for row, (rowax, df) in enumerate(zip(axs, dfs)):
            for ax, name_root in zip(rowax, module.c_name_roots):
                if module.logx == True:
                    ax.set_xscale('log', base=2)
                if module.logy == True:
                    ax.set_yscale('log', base=2)
                ax.set_aspect('equal', 'box')
                dif = df.loc[df.Time == t, name_root] - dfs[1].loc[df.Time == t, name_root]
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
                if name_root != module.c_name_roots[0]:
                    ax.set(yticks=[])
                ax.set_xlim(module.x_min, module.x_max)
                ax.set_ylim(module.y_min, module.y_max)
                ax.tick_params(axis='x', labelsize=fstick)
                ax.tick_params(axis='y', labelsize=fstick)
                #if name_root == 'a2Seen31':
                #    s = s + df.loc[df.Time == t, 'a2SeenSD31']
                #elif name_root == 'w33':
                #    s = s + df.loc[df.Time == t, 'wSD33']
                #elif name_root == 'w43':
                #    s = s + df.loc[df.Time == t, 'wSD43']
                #else:
                #    s = s + df.loc[df.Time == t, name_root + 'SD']
                s = s + df.loc[df.Time == t, name_root + 'SD']
                ax.scatter(x, y, c=color, edgecolor='0.400', alpha=0.2, s=s*module.bubble_size)
                if (row == 0):
                    ax.set(xticks=[])
                if name_root != module.c_name_roots[0]:
                    ax.set(yticks=[])
                ax.set_xlim(module.x_min, module.x_max)
                ax.set_ylim(module.y_min, module.y_max)
                ax.tick_params(axis='x', labelsize=fstick)
                ax.tick_params(axis='y', labelsize=fstick)

        plt.savefig(outfile)
        plt.close()

def get_data(dfs):

    dirs = ['', 'drift/']

    for d in dirs:
        df = pd.concat(map(pd.read_csv, glob.glob(os.path.join(d, '*.csv'))), ignore_index=True)
        if (sys.argv[1] == 'gsscattergrain') or (sys.argv[1] == 'gsscattergraini') or (sys.argv[1] == 'gsbars') or (sys.argv[1] == 'gsbarsone'):
            df = df[df.ChooseCost == 0.000061]
        if module.movie == False:
            lastt = df.Time.iat[-1]
            df = df[df.Time == lastt]
        dfs.append(df)

    return dfs

dfs = []
dfs = get_data(dfs)

if module.ftype == 'barsone':
    pr = barsonepr()
elif module.ftype == 'barsall':
    pr = barsallpr()
else:
    pr = scatterpr()

pr.prepare(dfs)

if module.movie == True:
    outfiles = []
    for t in dfs[0].Time.unique():
        print(f'Processing time {t}', end='\r')
        outfile = f'delete{t}.png'
        pr.chart(dfs, t)
        outfiles.append(outfile)
    giffile = f'{sys.argv[1]}.gif'
    with imageio.get_writer(giffile, mode='I') as writer:
        for outfile in outfiles:
            print(f'Adding {outfile} to movie', end='\r')
            image = imageio.imread(outfile)
            writer.append_data(image)
    for outfile in set(outfiles):
        os.remove(outfile)
else:
    pr.chart(dfs, dfs[0].Time.iat[-1])
