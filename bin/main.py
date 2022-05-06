#!/usr/bin/python

import glob
import imageio
import importlib.util
import math
import os
import sys
import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if len(sys.argv) < 2:
    print('You must add an argument with the python module name (for example, scatter)')
    exit(1)

start_time = time.perf_counter ()

sys.path.append(os.path.dirname(os.path.expanduser(sys.argv[1])))
module = __import__(sys.argv[1])

outfile = module.filename + '.png'

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

width = 12.0
height = 6.2
fslabel = 18 # Label font size
fstitle= 18 # Label font size
fstick = 14 # Label font size
red = 0.97
green = 0.97
blue = 0.97

lst_xy = [('ES', 'Substitutability of resources', True, None, None), 
            ('Given', 'Proportion of resource $\it{A}$\ngiven to partner', False, None, None),
            ('ChooseCost', 'Cost of comparing potential partners', True, None, None),
            ('MimicCost', 'Cost of comparing partner to self', True, None, None),
            ('DeathRate', 'Death rate', True, 0.005, 0.2),
            ('GroupSize', 'Number of potential partners', True, 48.0, 2.6),
            ('MimicGrainInit', 'Sensitivity for comparing partner to self', False, None, None),
            ('ChooseGrainInit', 'Sensitivity for comparing potential partners', False, None, None)]

lst_z = [('ChooseGrainmedian', 'Sensitivity for comparing\npotential partners', 600.0, None),
            ('ChooseGrain', 'Sensitivity for comparing\npotential partners', 600.0, 0.20),
            ('MimicGrainmedian', 'Sensitivity for comparing\npartner to self', 600.0, None),
            ('MimicGrain', 'Sensitivity for comparing\npartner to self', 600.0, 0.20),
            ('helpmedian', 'Help', 600.0, None),
            ('help', 'Help', 600.0, 0.20),
            ('wmedian', 'Fitness', 600.0, None),
            ('w', 'Fitness', 600.0, 0.20),
            ('chose_partner', 'Frequency of\nswitching to a new partner', 2000.0, None),
            ('changed_a2', 'Frequency of\nchanging help', 2000.0, None),
            ('helpBD', 'Fluctuation of help', 2000.0, None),
            ('wBD', 'Fluctuation of fitness', 2000.0, None)]

df_xy = pd.DataFrame(lst_xy, columns = ['xy', 'label', 'log', 'xymin', 'xymax'])
df_z = pd.DataFrame(lst_z, columns = ['z', 'title', 'bubble_size', 'ymax'])

x_name = module.x
x_label = df_xy.loc[df_xy.xy == x_name, 'label'].values[0]
if pd.isnull(df_xy.loc[df_xy.xy == x_name, 'xymin'].values[0]):
    x_min = None
else:
    x_min = df_xy.loc[df_xy.xy == x_name, 'xymin'].values[0]
if pd.isnull(df_xy.loc[df_xy.xy == x_name, 'xymax'].values[0]):
    x_max = None
else:
    x_max = df_xy.loc[df_xy.xy == x_name, 'xymax'].values[0]
x_log = df_xy.loc[df_xy.xy == x_name, 'log'].values[0]
y_name = module.y
y_label = df_xy.loc[df_xy.xy == y_name, 'label'].values[0]
if pd.isnull(df_xy.loc[df_xy.xy == y_name, 'xymin'].values[0]):
    y_min = None
else:
    y_min = df_xy.loc[df_xy.xy == y_name, 'xymin'].values[0]
if pd.isnull(df_xy.loc[df_xy.xy == y_name, 'xymax'].values[0]):
    y_max = None
else:
    y_max = df_xy.loc[df_xy.xy == y_name, 'xymax'].values[0]
y_log = df_xy.loc[df_xy.xy == y_name, 'log'].values[0]
z0_name = module.z0
z1_name = module.z1
c_name_roots = [z0_name, z1_name]
titles = [df_z.loc[df_z.z == z0_name, 'title'].values[0], df_z.loc[df_z.z == z1_name, 'title'].values[0]]

class barsallpr:

    def prepare(self, dfs):

        self.bincount = int(sum(map(lambda x: c_name_roots[0] in x, [*dfs[0]]))/2) - 2

        self.bh_maxs = [df_z.loc[df_z.z == z0_name, 'ymax'].values[0], df_z.loc[df_z.z == z1_name, 'ymax'].values[0]]

        self.c_name_suffixes = [x for x in range(self.bincount)]
        self.c_names = []
        for root in c_name_roots:
            rootlist = []
            for suffix in self.c_name_suffixes:
                rootlist.append(root + str(suffix)) 
            self.c_names.append(rootlist)

        self.inner_cols = dfs[0][x_name].unique()
        self.inner_rows = dfs[0][y_name].unique()
        self.inner_cols.sort()
        if y_name == 'GroupSize':
            self.inner_rows.sort()
        else:
            self.inner_rows[::-1].sort()

        self.bins = [(x+1)/self.bincount for x in range(self.bincount)]

        return self

    def chart(self, dfs, t):
 
        fig = plt.figure(figsize=(width, height-0.5))
        fig.supxlabel(x_label, fontsize=fslabel)
        fig.supylabel(t=y_label, x=0.003*width, fontsize=fslabel, ha='center')

        if module.movie == True:
            fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        outer_grid = fig.add_gridspec(nrows=1, ncols=len(c_name_roots), wspace=0.1)

        barwidth = 2.0/self.bincount

        for n, (c_name, name_root, title, bh_max) in enumerate(zip(self.c_names, c_name_roots, titles, self.bh_maxs)):
            inner_grid = outer_grid[n].subgridspec(nrows=len(self.inner_rows), ncols=len(self.inner_cols), wspace=0.0, hspace=0.0)
            axs = inner_grid.subplots()
            axs[0, int(len(self.inner_cols)/2)].set_title(title, fontsize=fslabel) # Prints the title of the middle column. Bad if there are even columns
            for row, (rowax, inner_row) in enumerate(zip(axs, self.inner_rows)): 
                for column, (ax, inner_col) in enumerate(zip(rowax, self.inner_cols)):
                    df = dfs[0]
                    median0 = df.loc[(df[x_name] == inner_col) & (df[y_name] == inner_row) & (df.Time == t), name_root + 'median'].values[0]
                    df = dfs[1]
                    median1 = df.loc[(df[x_name] == inner_col) & (df[y_name] == inner_row) & (df.Time == t), name_root + 'median'].values[0]
                    dif = median0 - median1
                    if (name_root == 'ChooseGrain') or (name_root == 'MimicGrain') or ('BD' in name_root):
                        dif = -dif
                    if dif < 0.0:
                        color=(red-0.95, green-0.95, blue-0.05)
                    else:
                        color=(red-0.95, green-0.05, blue-0.95)
                    df = dfs[0]
                    for b, name0, name1 in zip(self.bins[::2], c_name[::2], c_name[1::2]):
                        barheight=df.loc[(df[x_name] == inner_col) & (df[y_name] == inner_row) & (df.Time == t), name0] + df.loc[(df[x_name] == inner_col) & (df[y_name] == inner_row) & (df.Time == t), name1]
                        ax.bar(x=b, height=barheight, align='edge', color=color, linewidth=0, width=barwidth)
                        ax.set(xticks=[], yticks=[], ylim=[0, bh_max*2])
                    df = dfs[1]
                    for b, name0, name1 in zip(self.bins[::2], c_name[::2], c_name[1::2]):
                        barheight=df.loc[(df[x_name] == inner_col) & (df[y_name] == inner_row) & (df.Time == t), name0] + df.loc[(df[x_name] == inner_col) & (df[y_name] == inner_row) & (df.Time == t), name1]
                        ax.bar(x=b, height=barheight, align='edge', color=(red-0.15, green-0.15, blue-0.15), linewidth=0, width=barwidth, alpha=0.9)
                        ax.set(xticks=[], yticks=[], ylim=[0, bh_max*2])
                    if (n == 0) & (column == 0):
                        if y_log == True:
                            y = '$2^{{{}}}$'.format(round(math.log(inner_row, 2)))
                        else:
                            y = inner_row
                        ax.set_ylabel(y, rotation='horizontal', horizontalalignment='right', verticalalignment='center')
                    if row == len(self.inner_rows) - 1:
                        if x_log == True:
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

        if x_log == True:
            self.x_value = float(str("{:.6f}".format(pow(2, int(module.x_value)))))
        else:
            self.x_value = module.x_value
        if y_log == True:
            self.y_value = float(str("{:.6f}".format(pow(2, int(module.y_value)))))
        else:
            self.y_value = module.y_value

        self.bh_maxs = [df_z.loc[df_z.z == z0_name, 'ymax'].values[0], df_z.loc[df_z.z == z1_name, 'ymax'].values[0]]

        self.bincount = int(sum(map(lambda x: c_name_roots[0] in x, [*dfs[0]]))/2) - 2

        self.c_name_suffixes = [x for x in range(self.bincount)]
        self.c_names = []
        self.c_names_sd = []
        for root in c_name_roots:
            rootlist = []
            rootlist_sd = []
            for suffix in self.c_name_suffixes:
                rootlist.append(root + str(suffix)) 
                rootlist_sd.append(root + 'SD' + str(suffix)) 
            self.c_names.append(rootlist)
            self.c_names_sd.append(rootlist_sd)

        self.bins = [(x+1)/self.bincount for x in range(self.bincount)]

        return self

    def chart(self, dfs, t):

        fig, axs = plt.subplots(nrows=1, ncols=len(c_name_roots), figsize=(width, height), constrained_layout=False)

        fig.supylabel('\nFrequency', fontsize=fslabel, ha='center')

        if module.movie == True:
            fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        barwidth = -1.0/self.bincount

        for ax, title in zip(axs, titles):
            ax.set_xlabel(title, fontsize=fslabel)

        for ax, name_root, c_name, c_name_sd, bh_max in zip(axs, c_name_roots, self.c_names, self.c_names_sd, self.bh_maxs):
            df = dfs[0]
            median0 = df.loc[(df[x_name] == self.x_value) & (df[y_name] == self.y_value) & (df.Time == t), name_root + 'median'].values[0]
            df = dfs[1]
            median1 = df.loc[(df[x_name] == self.x_value) & (df[y_name] == self.y_value) & (df.Time == t), name_root + 'median'].values[0]
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
                barheight = df.loc[(df[x_name] == self.x_value) & (df[y_name] == self.y_value) & (df.Time == t), name]
                barheightsd = df.loc[(df[x_name] == self.x_value) & (df[y_name] == self.y_value) & (df.Time == t), namesd]
                ax.bar(x=b, height=barheight, align='edge', color=color, linewidth=0, width=barwidth, alpha=alpha)
                ax.bar(x=b, height=barheightsd, align='edge', color=colorsd, linewidth=0, width=barwidth, bottom=barheight, alpha=alpha)
            ax.set(ylim=(0, bh_max), yticks=(0, bh_max), yticklabels=(0, bh_max))
            ax.tick_params(axis='x', labelsize=fstick)
            ax.tick_params(axis='y', labelsize=fstick)
            if name_root != c_name_roots[0]:
                ax.set(yticks=[])
            ax.set_box_aspect(1)
            color=(red-0.15, green-0.15, blue-0.15)
            colorsd=(red-0.10, green-0.10, blue-0.10)
            alpha=0.9
            df = dfs[1]
            for b, name, namesd in zip(self.bins, c_name, c_name_sd):
                barheight = df.loc[(df[x_name] == self.x_value) & (df[y_name] == self.y_value) & (df.Time == t), name]
                barheightsd = df.loc[(df[x_name] == self.x_value) & (df[y_name] == self.y_value) & (df.Time == t), namesd]
                ax.bar(x=b, height=barheight, align='edge', color=color, linewidth=0, width=barwidth, alpha=alpha)
                ax.bar(x=b, height=barheightsd, align='edge', color=colorsd, linewidth=0, width=barwidth, bottom=barheight, alpha=alpha)
            ax.set(ylim=(0, bh_max), yticks=(0, bh_max), yticklabels=(0, bh_max))
            ax.tick_params(axis='x', labelsize=fstick)
            ax.tick_params(axis='y', labelsize=fstick)
            if name_root != c_name_roots[0]:
                ax.set(yticks=[])
            ax.set_box_aspect(1)

        plt.savefig(outfile, dpi=100)
        plt.close()

class scatterpr:

    def prepare(self, dfs):
        self.suffixes = ('SD', '')
        self.suffixalphas = (0.2, 1.0)
        self.suffixecs = ('0.300', '0.000')
        self.bubble_sizes = [df_z.loc[df_z.z == module.z0, 'bubble_size'].values[0], df_z.loc[df_z.z == module.z1, 'bubble_size'].values[0]]

    def chart(self, dfs, t):

        fig, axs = plt.subplots(nrows=1, ncols=len(c_name_roots), figsize=(width, height), constrained_layout=False)

        fig.supxlabel(x_label, fontsize=fslabel)
        fig.supylabel(y_label, fontsize=fslabel, ha='center')

        if module.movie == True:
            fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        for ax, name_root, title, bubble_size in zip(axs, c_name_roots, titles, self.bubble_sizes):
            ax.set_title(title, fontsize=fstitle)
            ax.tick_params(axis='x', labelsize=fstick)
            ax.tick_params(axis='y', labelsize=fstick)
            if x_log == True:
                ax.set_xscale('log', base=2)
            if y_log == True:
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
            for suffix, suffixec, suffixalpha in zip(self.suffixes, self.suffixecs, self.suffixalphas):
                df = dfs[module.drift]
                x = df.loc[df.Time == t, x_name]
                y = df.loc[df.Time == t, y_name]
                s = df.loc[df.Time == t, name_root]
                if suffix == 'SD':
                    s = s + df.loc[df.Time == t, name_root + suffix]
                if module.drift == True:
                    ax.scatter(x, y, color='0.700', edgecolor='0.700', alpha=suffixalpha, s=s*bubble_size)
                else:
                    ax.scatter(x, y, c=color, ec=color, alpha=suffixalpha, s=s*bubble_size)
                if name_root != c_name_roots[0]:
                    ax.set(yticks=[])
                ax.set_xlim(x_min, x_max)
                ax.set_ylim(y_min, y_max)
                ax.set_box_aspect(1)

        plt.savefig(outfile, transparent=False)
        plt.close()

def get_data(dfs):

    dirs = [module.treatment, module.control]

    for d in dirs:
        df = pd.concat(map(pd.read_csv, glob.glob(os.path.join(d, '*.csv'))), ignore_index=True)
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
    giffile = module.filename + '.gif'
    with imageio.get_writer(giffile, mode='I') as writer:
        for outfile in outfiles:
            print(f'Adding {outfile} to movie', end='\r')
            image = imageio.imread(outfile)
            writer.append_data(image)
    for outfile in set(outfiles):
        os.remove(outfile)
else:
    pr.chart(dfs, dfs[0].Time.iat[-1])

end_time = time.perf_counter ()
print(end_time - start_time, "seconds")
