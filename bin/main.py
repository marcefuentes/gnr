#!/usr/bin/python

import glob
import importlib.util
import math
import os
import sys
import time
import imageio.v3 as iio
import matplotlib.pyplot as plt
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
fstitle= 18 # Title font size
fstick = 14 # Tick font size
red = 0.97
green = 0.97
blue = 0.97

lst_xy = [('ES', 'Substitutability of resources', True, pow(2, -5.5), None), 
            ('Given', 'Proportion of resource $\it{A}$\ngiven to partner', False, None, None),
            ('ChooseCost', 'Cost of comparing potential partners', True, None, None),
            ('MimicCost', 'Cost of comparing partner to self', True, None, None),
            ('DeathRate', 'Death rate', True, 0.005, 0.2),
            ('GroupSize', 'Number of potential partners', True, 48.0, 2.6),
            ('MimicGrainInit', 'Sensitivity for comparing partner to self', False, None, None),
            ('ChooseGrainInit', 'Sensitivity for comparing potential partners', False, None, None),
            ('N', 'Population size', True, None, None)]

lst_z = [('ChooseGrainmedian', 'Sensitivity for comparing\npotential partners', 600.0, None),
            ('ChooseGrain', 'Sensitivity for comparing\npotential partners', 600.0, 0.20),
            ('MimicGrainmedian', 'Sensitivity for comparing\npartner to self', 600.0, None),
            ('MimicGrain', 'Sensitivity for comparing\npartner to self', 600.0, 0.20),
            ('helpmedian', 'Help', 600.0, None),
            ('help', 'Help', 600.0, 0.20),
            ('a2Seenmedian', 'a2', 600.0, None),
            ('a2Seen', 'a2', 600.0, 0.20),
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
z_names = [module.z0, module.z1]

titles = []
for z_name in z_names: 
    titles.append(df_z.loc[df_z.z == z_name, 'title'].values[0])

class barsallpr:

    def prepare(self, dfs):

        bincount = int(sum(map(lambda x: z_names[0] in x, [*dfs[0]]))/2) - 2

        self.z_namebins_lists = []
        self.bh_maxs = []
        for z_name in z_names:
            self.z_namebins_lists.append([z_name + str(x) for x in range(bincount)])
            self.bh_maxs.append(df_z.loc[df_z.z == z_name, 'ymax'].values[0])

        self.inner_cols = dfs[0][x_name].unique()
        self.inner_rows = dfs[0][y_name].unique()
        self.inner_cols.sort()
        if y_name == 'GroupSize':
            self.inner_rows.sort()
        else:
            self.inner_rows[::-1].sort()

        self.bins = [(x+1)/bincount for x in range(bincount)]
        self.barwidth = 2.0/bincount

        self.color_blue = [red-0.95, green-0.95, blue-0.05]
        self.color_green = [red-0.95, green-0.05, blue-0.95]
        color_gray = [red-0.15, green-0.15, blue-0.15]
        colorsd_gray = [red-0.10, green-0.10, blue-0.10]

        self.colors = [self.color_blue, color_gray]
        self.alphas = [1.0, 0.9]

        return self

    def chart(self, dfts):
 
        fig = plt.figure(figsize=(width, height-0.5))
        fig.supxlabel(x_label, fontsize=fslabel)
        fig.supylabel(t=y_label, x=0.003*width, fontsize=fslabel, ha='center')

        if module.movie == True:
            fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        outer_grid = fig.add_gridspec(nrows=1, ncols=len(z_names), wspace=0.1)

        for n, (z_namebins, z_name, title, bh_max) in enumerate(zip(self.z_namebins_lists, z_names, titles, self.bh_maxs)):
            inner_grid = outer_grid[n].subgridspec(nrows=len(self.inner_rows), ncols=len(self.inner_cols), wspace=0.0, hspace=0.0)
            axs = inner_grid.subplots()
            axs[0, int(len(self.inner_cols)/2)].set_title(title, fontsize=fslabel) # Prints the title of the middle column. Bad if there are even columns
            for row, (rowax, inner_row) in enumerate(zip(axs, self.inner_rows)): 
                for column, (ax, inner_col) in enumerate(zip(rowax, self.inner_cols)):
                    medians = []
                    for df in dfts:
                        medians.append(df.loc[(df[x_name] == inner_col) & (df[y_name] == inner_row), z_name + 'median'].values[0])
                    dif = medians[0] - medians[1]
                    if ('Grain' in z_name) or ('BD' in z_name):
                        dif = -dif
                    if dif > 0.0:
                        self.colors[0] = self.color_green
                    else:
                        self.colors[0] = self.color_blue

                    for df, color, alpha in zip(dfts, self.colors, self.alphas):
                        for b, name0, name1 in zip(self.bins[::2], z_namebins[::2], z_namebins[1::2]):
                            barheight=df.loc[(df[x_name] == inner_col) & (df[y_name] == inner_row), name0] + df.loc[(df[x_name] == inner_col) & (df[y_name] == inner_row), name1]
                            ax.bar(x=b, height=barheight, align='edge', color=color, linewidth=0, width=self.barwidth, alpha=alpha)
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

        bincount = int(sum(map(lambda x: z_names[0] in x, [*dfs[0]]))/2) - 2

        self.z_namebins_lists = []
        self.z_namesdbins_lists = []
        self.bh_maxs = []
        self.binslists = []
        self.barwidths = []
        for z_name in z_names:
            self.z_namebins_lists.append([z_name + str(x) for x in range(bincount)])
            self.z_namesdbins_lists.append([z_name + 'SD' + str(x) for x in range(bincount)])
            self.bh_maxs.append(df_z.loc[df_z.z == z_name, 'ymax'].values[0])
            if z_name == 'w':
                mmax = 2.0 # For a1Max = a2Max = 1.0 and R1 = R2 = 2.0.
            else:
                mmax = 1.0
            self.binslists.append([(x+1)*mmax/bincount for x in range(bincount)])
            self.barwidths.append(-mmax/bincount)

        self.color_blue = [red-0.95, green-0.95, blue-0.05]
        self.colorsd_blue = [red-0.30, green-0.30, blue-0.05]
        self.color_green = [red-0.95, green-0.05, blue-0.95]
        self.colorsd_green = [red-0.30, green-0.05, blue-0.30]
        color_gray = [red-0.15, green-0.15, blue-0.15]
        colorsd_gray = [red-0.10, green-0.10, blue-0.10]

        self.colors = [self.color_blue, color_gray]
        self.colorsds= [self.colorsd_blue, colorsd_gray]
        self.alphas = [1.0, 0.8]

        return self

    def chart(self, dfts):

        fig, axs = plt.subplots(nrows=1, ncols=len(z_names), figsize=(width, height), constrained_layout=False)

        fig.supylabel('\nFrequency', fontsize=fslabel, ha='center')

        if module.movie == True:
            fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        for ax, title in zip(axs, titles):
            ax.set_xlabel(title, fontsize=fslabel)

        for ax, z_name, z_namebins, z_namesdbins, bins, barwidth, bh_max in zip(axs, z_names, self.z_namebins_lists, self.z_namesdbins_lists, self.binslists, self.barwidths, self.bh_maxs):
            medians = []
            for df in dfts:
                medians.append(df.loc[(df[x_name] == self.x_value) & (df[y_name] == self.y_value), z_name + 'median'].values[0])
            dif = medians[0] - medians[1]
            if ('Grain' in z_name) or ('BD' in z_name):
                dif = -dif
            if dif > 0.0:
                self.colors[0] = self.color_green
                self.colorsds[0] = self.colorsd_green
            else:
                self.colors[0] = self.color_blue
                self.colorsds[0] = self.colorsd_blue
            for df, color, colorsd, alpha in zip(dfts, self.colors, self.colorsds, self.alphas):
                for b, z_namebin, z_namesdbin in zip(bins, z_namebins, z_namesdbins):
                    barheight = df.loc[(df[x_name] == self.x_value) & (df[y_name] == self.y_value), z_namebin]
                    barheightsd = df.loc[(df[x_name] == self.x_value) & (df[y_name] == self.y_value), z_namesdbin]
                    ax.bar(x=b, height=barheight, align='edge', color=color, linewidth=0, width=barwidth, alpha=alpha)
                    ax.bar(x=b, height=barheightsd, align='edge', color=colorsd, linewidth=0, width=barwidth, bottom=barheight, alpha=alpha)
                ax.set(ylim=(0, bh_max), yticks=(0, bh_max), yticklabels=(0, bh_max))
                ax.tick_params(axis='x', labelsize=fstick)
                ax.tick_params(axis='y', labelsize=fstick)
                if z_name != z_names[0]:
                    ax.set(yticks=[])
                ax.set_box_aspect(1)

        plt.savefig(outfile, dpi=100)
        plt.close()

class scatterpr:

    def prepare(self, dfs):
        self.suffixes = ('SD', '')
        self.suffixalphas = (0.2, 1.0)
        self.suffixecs = ('0.300', '0.000')
        self.bubble_sizes = []
        for z_name in z_names:
            self.bubble_sizes.append(df_z.loc[df_z.z == z_name, 'bubble_size'].values[0])

    def chart(self, dfts):

        fig, axs = plt.subplots(nrows=1, ncols=len(z_names), figsize=(width, height), constrained_layout=False)

        fig.supxlabel(t=x_label, y=0.02, fontsize=fslabel)
        fig.supylabel(t=y_label, x=0.04, fontsize=fslabel, ha='center')

        if module.movie == True:
            fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        for dft in dfts:
            dft.sort_values(by=[x_name, y_name], inplace=True)

        for ax, z_name, title, bubble_size in zip(axs, z_names, titles, self.bubble_sizes):
            ax.set_title(title, pad=10.0, fontsize=fstitle)
            ax.tick_params(axis='x', labelsize=fstick)
            ax.tick_params(axis='y', labelsize=fstick)
            if x_log == True:
                ax.set_xscale('log', base=2)
            if y_log == True:
                ax.set_yscale('log', base=2)
            dif = dfts[1][z_name] - dfts[0][z_name]
            if (z_name == 'ChooseGrainmedian') or (z_name == 'MimicGrainmedian') or ('BD' in z_name):
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
                df = dfts[module.drift]
                x = df[x_name]
                y = df[y_name]
                s = df[z_name]
                if suffix == 'SD':
                    s = s + df[z_name + suffix]
                if module.drift == True:
                    ax.scatter(x, y, color='0.700', edgecolor='0.700', alpha=suffixalpha, s=s*bubble_size)
                else:
                    ax.scatter(x, y, c=color, ec=color, alpha=suffixalpha, s=s*bubble_size)
                if z_name != z_names[0]:
                    ax.axes.yaxis.set_ticklabels([])
                ax.set_xlim(x_min, x_max)
                ax.set_ylim(y_min, y_max)
                ax.set_box_aspect(1)

        plt.savefig(outfile, transparent=False)
        plt.close()

def get_data(dfs):

    dirs = [module.treatment, module.control]

    for d in dirs:
        df = pd.concat(map(pd.read_csv, glob.glob(os.path.join(d, '*.csv'))), ignore_index=True)
        dfs.append(df)

    if module.control == 'optimal':
        dfs[1]['helpmedian'] = dfs[1]['a2Seenmedian']*2.0*dfs[0]['Given']
    if module.treatment == 'optimal':
        dfs[0]['helpmedian'] = dfs[0]['a2Seenmedian']*2.0*dfs[1]['Given']

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
        print(f'Processing step {t}', end='\r')
        outfile = f'delete{t}.png'
        dfts = []
        for df in dfs:
            dfts.append(df[df.Time == t].copy())
        pr.chart(dfts)
        outfiles.append(outfile)
    frames = []
    for outfile in outfiles:
        frames.append(iio.imread(outfile))
    giffile = module.filename + '.gif'
    iio.imwrite(giffile, frames)
    [os.remove(outfile) for outfile in set(outfiles)]
else:
    t = dfs[0].Time.iat[-1]
    dfts = []
    for df in dfs:
        dfts.append(df[df.Time == t].copy())
    pr.chart(dfts)

end_time = time.perf_counter ()
print("\nTime elapsed: %.2f seconds" % (end_time - start_time))
