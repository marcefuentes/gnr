#! /usr/bin/env python

import glob
import math
import os
import sys
import time
import imageio.v3 as iio
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

width = 6.0*len(module.z)
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

dfxy = pd.DataFrame(lst_xy, columns = ['xy', 'label', 'log', 'xymin', 'xymax'])
dfz = pd.DataFrame(lst_z, columns = ['z', 'title', 'bubble_size', 'ymax'])

xname = module.x
yname = module.y
x_label = dfxy.loc[dfxy.xy == xname, 'label'].values[0]
y_label = dfxy.loc[dfxy.xy == yname, 'label'].values[0]
x_log = dfxy.loc[dfxy.xy == xname, 'log'].values[0]
y_log = dfxy.loc[dfxy.xy == yname, 'log'].values[0]
x_min = None if pd.isnull(dfxy.loc[dfxy.xy == xname, 'xymin'].values[0]) else dfxy.loc[dfxy.xy == xname, 'xymin'].values[0]
y_min = None if pd.isnull(dfxy.loc[dfxy.xy == yname, 'xymin'].values[0]) else dfxy.loc[dfxy.xy == yname, 'xymin'].values[0]
x_max = None if pd.isnull(dfxy.loc[dfxy.xy == xname, 'xymax'].values[0]) else dfxy.loc[dfxy.xy == xname, 'xymax'].values[0]
y_max = None if pd.isnull(dfxy.loc[dfxy.xy == yname, 'xymax'].values[0]) else dfxy.loc[dfxy.xy == yname, 'xymax'].values[0]

class BarsAll:

    def prepare(self, dfs):

        bincount = int(sum(map(lambda x: module.z0['name'] in x, [*dfs[module.dirs[0]]]))/2) - 2

        for zdict in zdicts:
            zdict['namebins_list'] = [zdict['name'] + str(x) for x in range(bincount)]
            zdict['bh_max'] = dfz.loc[dfz.z == zdict['name'], 'ymax'].values[0]

        self.innercols = dfs[module.dirs[0]][xname].unique()
        self.innerrows = dfs[module.dirs[0]][yname].unique()
        self.innercols.sort()
        self.innerrows.sort() if yname == 'GroupSize' else self.innerrows[::-1].sort()

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

        if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        outer_grid = fig.add_gridspec(nrows=1, ncols=len(zdicts), wspace=0.1)

        for n, zdict in enumerate(zdicts):
            innergrid = outer_grid[n].subgridspec(nrows=len(self.innerrows), ncols=len(self.innercols), wspace=0.0, hspace=0.0)
            axs = innergrid.subplots()
            axs[0, int(len(self.innercols)/2)].set_title(zdict['title'], fontsize=fslabel) # Prints the title of the middle column. Bad if there are even columns
            for row, (rowax, innerrow) in enumerate(zip(axs, self.innerrows)): 
                for column, (ax, innercol) in enumerate(zip(rowax, self.innercols)):
                    ds = [dfts[zdict['control']], dfts[zdict['treatment']]]
                    medians = []
                    [medians.append(d.loc[(d[xname] == innercol) & (d[yname] == innerrow), zdict['name'] + 'median'].values[0]) for d in ds]
                    dif = medians[1] - medians[0]
                    if ('Grain' in zdict['name']) or ('BD' in zdict['name']): dif = -dif
                    self.colors[0] = self.color_green if dif > 0.0 else self.color_blue
                    for d, color, alpha in zip(ds, self.colors, self.alphas):
                        for b, name0, name1 in zip(self.bins[::2], zdict['namebins_list'][::2], zdict['namebins_list'][1::2]):
                            barheight=d.loc[(d[xname] == innercol) & (d[yname] == innerrow), name0] + d.loc[(d[xname] == innercol) & (d[yname] == innerrow), name1]
                            ax.bar(x=b, height=barheight, align='edge', color=color, linewidth=0, width=self.barwidth, alpha=alpha)
                            ax.set(xticks=[], yticks=[], ylim=[0, zdict['bh_max']*2])
                    if (n == 0) & (column == 0):
                        y = '$2^{{{}}}$'.format(round(math.log(innerrow, 2))) if y_log else innerrow
                        ax.set_ylabel(y, rotation='horizontal', horizontalalignment='right', verticalalignment='center')
                    if row == len(self.innerrows) - 1:
                        x = '$2^{{{}}}$'.format(round(math.log(innercol, 2))) if x_log else innercol
                        ax.set_xlabel(x)

        plt.savefig(outfile, dpi=100)
        plt.close()

class BarsOne:

    def prepare(self, dfs):

        self.x_value = float(str("{:.6f}".format(pow(2, int(module.x_value))))) if x_log else module.x_value
        self.y_value = float(str("{:.6f}".format(pow(2, int(module.y_value))))) if y_log else module.y_value

        bincount = int(sum(map(lambda x: module.z0['name'] in x, [*dfs[module.dirs[0]]]))/2) - 2

        for zdict in zdicts:
            zdict['namebins_list'] = [zdict['name'] + str(x) for x in range(bincount)]
            zdict['namesdbins_list'] = [zdict['name'] + 'SD' + str(x) for x in range(bincount)]
            zdict['bh_max'] = dfz.loc[dfz.z == zdict['name'], 'ymax'].values[0]
            mmax = 2.0 if zdict['name'] == 'w' else 1.0 # For a1Max = a2Max = 1.0 and R1 = R2 = 2.0.
            zdict['binslist'] = [(x+1)*mmax/bincount for x in range(bincount)]
            zdict['barwidth'] = -mmax/bincount

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

        fig, axs = plt.subplots(nrows=1, ncols=len(zdicts), figsize=(width, height), sharey=True, constrained_layout=False, squeeze=False)

        fig.supylabel('\nFrequency', fontsize=fslabel, ha='center')

        if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        for ax, zdict in zip(axs.reshape(-1), zdicts):
            ax.set_xlabel(zdict['title'], fontsize=fslabel)
            ds = [dfts[zdict['control']], dfts[zdict['treatment']]]
            medians = []
            [medians.append(d.loc[(d[xname] == self.x_value) & (d[yname] == self.y_value), zdict['name'] + 'median'].values[0]) for d in ds]
            dif = medians[1] - medians[0]
            if ('Grain' in zdict['name']) or ('BD' in zdict['name']): dif = -dif
            self.colors[0] = self.color_green if dif > 0.0 else self.color_blue
            self.colorsds[0] = self.color_green if dif > 0.0 else self.colorsd_blue
            for d, color, colorsd, alpha in zip(ds, self.colors, self.colorsds, self.alphas):
                for b, namebin, namesdbin in zip(zdict['binslist'], zdict['namebins_list'], zdict['namesdbins_list']):
                    barheight = d.loc[(d[xname] == self.x_value) & (d[yname] == self.y_value), namebin]
                    barheightsd = d.loc[(d[xname] == self.x_value) & (d[yname] == self.y_value), namesdbin]
                    ax.bar(x=b, height=barheight, align='edge', color=color, linewidth=0, width=zdict['barwidth'], alpha=alpha)
                    ax.bar(x=b, height=barheightsd, align='edge', color=colorsd, linewidth=0, width=zdict['barwidth'], bottom=barheight, alpha=alpha)
                ax.set(ylim=(0, zdict['bh_max']), yticks=(0, zdict['bh_max']), yticklabels=(0, zdict['bh_max']))
                ax.tick_params(axis='x', labelsize=fstick)
                ax.tick_params(axis='y', labelsize=fstick)
                ax.set_box_aspect(1)

        plt.savefig(outfile, dpi=100)
        plt.close()

class Scatter:

    def prepare(self, dfs):
        self.suffixes = ('SD', '')
        self.suffixalphas = (0.2, 1.0)
        self.suffixecs = ('0.300', '0.000')
        for zdict in zdicts:
            zdict['bubble_size'] = dfz.loc[dfz.z == zdict['name'], 'bubble_size'].values[0]

    def chart(self, dfts):

        fig, axs = plt.subplots(nrows=1, ncols=len(zdicts), figsize=(width, height), sharey=True, constrained_layout=False, squeeze=False)
        fig.supxlabel(t=x_label, y=0.02, fontsize=fslabel)
        fig.supylabel(t=y_label, x=0.04, fontsize=fslabel, ha='center')

        if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        for d in module.dirs:
            dfts[d].sort_values(by=[xname, yname], inplace=True)

        for ax, zdict in zip(axs.reshape(-1), zdicts):
            ax.set_title(zdict['title'], pad=10.0, fontsize=fstitle)
            ax.tick_params(axis='x', labelsize=fstick)
            ax.tick_params(axis='y', labelsize=fstick)
            if x_log: ax.set_xscale('log', base=2)
            if y_log: ax.set_yscale('log', base=2)
            df = dfts[zdict['treatment']]
            x = df[xname]
            y = df[yname]
            s = df[zdict['name']]
            dif = dfts[zdict['control']][zdict['name']] - s
            if (zdict['name'] == 'ChooseGrainmedian') or (zdict['name'] == 'MimicGrainmedian') or ('BD' in zdict['name']): dif = -dif
            color = []
            for i in dif:
                if i < 0.0:
                    if i < -red: i = -red
                    color.append((red + i, green + i/2.0, blue + i))
                else:
                    if i > red: i = red
                    color.append((red - i, green - i, blue - i/2.0))
            for suffix, suffixec, suffixalpha in zip(self.suffixes, self.suffixecs, self.suffixalphas):
                if suffix == 'SD': s = s + df[zdict['name'] + suffix]
                if zdict['treatment'] == 'none':
                    ax.scatter(x, y, color='0.700', edgecolor='0.700', alpha=suffixalpha, s=s*zdict['bubble_size'])
                else:
                    ax.scatter(x, y, c=color, ec=color, alpha=suffixalpha, s=s*zdict['bubble_size'])
                ax.set_xlim(x_min, x_max)
                ax.set_ylim(y_min, y_max)
                ax.set_box_aspect(1)

        plt.savefig(outfile, transparent=False)
        plt.close()

def create_figure(t):
    dfts = {}
    for d in module.dirs:
        dfts[d] = dfs[d].loc[dfs[d]['Time'] == t].copy()
    pr.chart(dfts)

zdicts = module.z

for zdict in zdicts:
    zdict['title'] = dfz.loc[dfz.z == zdict['name'], 'title'].values[0]

dfs = {}
for d in module.dirs:
    dfs[d] = pd.concat(map(pd.read_csv, glob.glob(os.path.join(d, '*.csv'))), ignore_index=True)
    if d == 'optimal': dfs[d]['helpmedian'] = dfs[d]['a2Seenmedian']*2.0*dfs[d]['Given']

if module.ftype == 'barsone':
    pr = BarsOne()
elif module.ftype == 'barsall':
    pr = BarsAll()
elif module.ftype == 'scatter':
    pr = Scatter()
else:
    print('No such ftype')
    exit(1)

pr.prepare(dfs)

if module.movie:
    frames = []
    for t in dfs[module.dirs[0]].Time.unique():
        print(f'Processing step {t}', end='\r')
        outfile = f'delete{t}.png'
        create_figure(t)
        frames.append(iio.imread(outfile))
        os.remove(outfile)
    giffile = module.filename + '.gif'
    iio.imwrite(giffile, frames)
else:
    t = dfs[module.dirs[0]].Time.iat[-1]
    outfile = module.filename + '.png'
    create_figure(t)

end_time = time.perf_counter ()
print("\nTime elapsed: %.2f seconds" % (end_time - start_time))
