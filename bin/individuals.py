#! /usr/bin/env python

import glob
import math
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

height = 6.0*len(module.z)
width = 6.0*len(module.z[0])
fslabel = 24 # Label font size
fstitle= 24 # Title font size
fstick = 16 # Tick font size
red = 0.97
green = 0.97
blue = 0.97

lst_xy = [('ES', 'Substitutability of resource $\it{A}$', True, pow(2, -5.5), None), 
            ('Given', 'Partner\'s share of resource $\it{A}$', False, None, None),
            ('ChooseCost', 'Cost of comparing potential partners', True, None, None),
            ('MimicCost', 'Cost of comparing partner to self', True, None, None),
            ('DeathRate', 'Death rate', True, 0.005, 0.2),
            ('GroupSize', 'Number of potential partners', True, 48.0, 2.6),
            ('MimicGrainInit', 'Sensitivity for comparing partner to self', False, None, None),
            ('ChooseGrainInit', 'Sensitivity for comparing potential partners', False, None, None),
            ('N', 'Population size', True, None, None)]

lst_z = [('ChooseGrainmedian', 'Sensitivity for comparing\npotential partners', 600.0, 0.0, 1.0),
            ('ChooseGrain', 'Sensitivity for comparing\npotential partners', 600.0, 0.0, 1.0),
            ('MimicGrainmedian', 'Sensitivity for comparing\npartner to self', 600.0, 0.0, 1.0),
            ('MimicGrain', 'Sensitivity for comparing\npartner to self', 600.0, 0.0, 1.0),
            ('helpmedian', 'Help', 600.0*1.87, 0.0, 1.87),
            ('help', 'Help', 600.0*1.87, 0.0, 1.87),
            ('a2Seenmedian', '$\it{a}$', 600.0, 0.0, 1.0),
            ('a2Seen', '$\it{a}$', 600.0, 0.0, 1.0),
            ('a2Defaultmedian', 'Default $\it{a}$', 600.0, 0.0, 1.0),
            ('a2Default', 'Default $\it{a}$', 600.0, 0.0, 1.0),
            ('wmedian', 'Fitness', 600.0, 0.0, 1.87),
            ('w', 'Fitness', 600.0, 0.0, 1.87),
            ('chose_partner', 'Frequency of\nswitching to a new partner', 2000.0, -0.1, 1.1),
            ('changed_a2', 'Frequency of\nchanging $\it{a}$', 2000.0, -0.1, 1.1),
            ('helpBD', 'Fluctuation of help', 2000.0, 0.0, 1.0),
            ('wBD', 'Fluctuation of fitness', 2000.0, 0.0, 1.0)]

dfxy = pd.DataFrame(lst_xy, columns = ['xy', 'label', 'log', 'xymin', 'xymax'])
dfz = pd.DataFrame(lst_z, columns = ['z', 'title', 'bubble_size', 'zmin', 'zmax'])

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

class ScatterAll:

    def prepare(self, dfs):

        for zdicts in zdictss:
            for zdict in zdicts:
                zdict['xmax'] = 2.0 if zdict['x'] == 'w' else 1.0 # For a1Max = a2Max = 1.0 and R1 = R2 = 2.0.
                zdict['ymax'] = 2.0 if zdict['y'] == 'w' else 1.0 # For a1Max = a2Max = 1.0 and R1 = R2 = 2.0.
                zdict['title'] = 'x = ' + dfz.loc[dfz.z == zdict['x'], 'title'].values[0] + '\ny = ' + dfz.loc[dfz.z == zdict['y'], 'title'].values[0]
        self.innercols = dfs[module.dirs[0]][xname].unique()
        self.innerrows = dfs[module.dirs[0]][yname].unique()
        self.innercols.sort()
        self.innerrows.sort() if yname == 'GroupSize' else self.innerrows[::-1].sort()

#        self.color = [red-0.15, green-0.15, blue-0.15]

        return self

    def chart(self, dfts):
 
        fig = plt.figure(figsize=(width, height))
        fig.supxlabel(x_label, fontsize=fslabel)
        fig.supylabel(t=y_label, x=0.003*width, fontsize=fslabel, ha='center')

        if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        outer_grid = fig.add_gridspec(nrows=len(zdictss), ncols=len(zdictss[0]), wspace=0.1)

        for nr, zdicts in enumerate(zdictss):
            for nc, zdict in enumerate(zdicts):
                dft = dfts[zdict['treatment']]
                innergrid = outer_grid[nr, nc].subgridspec(nrows=len(self.innerrows), ncols=len(self.innercols), wspace=0.0, hspace=0.0)
                axs = innergrid.subplots()
                if nr == 0:    
                    axs[0, int(len(self.innercols)/2)].set_title(zdict['title'], fontsize=fslabel) # Prints the title of the middle column. Bad if there are even columns
                for row, (rowax, innerrow) in enumerate(zip(axs, self.innerrows)): 
                    for column, (ax, innercol) in enumerate(zip(rowax, self.innercols)):
                        x = dft.loc[(dft[xname] == innercol) & (dft[yname] == innerrow), zdict['x']]
                        y = dft.loc[(dft[xname] == innercol) & (dft[yname] == innerrow), zdict['y']]
                        ax.scatter(x, y, alpha=0.1, s=0.001)
                        ax.set_xlim(0.0, zdict['xmax'])
                        ax.set_ylim(0.0, zdict['ymax'])
                        ax.set(xticks=[], yticks=[])
                        if (nc == 0) & (column == 0):
                            y = '$2^{{{}}}$'.format(round(math.log(innerrow, 2))) if y_log else innerrow
                            ax.set_ylabel(y, rotation='horizontal', horizontalalignment='right', verticalalignment='center')
                        if row == len(self.innerrows) - 1:
                            x = '$2^{{{}}}$'.format(round(math.log(innercol, 2))) if x_log else innercol
                            ax.set_xlabel(x)

        plt.savefig(outfile, dpi=100)
        plt.close()

class ScatterOne:

    def prepare(self, dfs):

        self.x_value = float(str("{:.6f}".format(pow(2, int(module.x_value))))) if x_log else module.x_value
        self.y_value = float(str("{:.6f}".format(pow(2, int(module.y_value))))) if y_log else module.y_value

        for zdicts in zdictss:
            for zdict in zdicts:
                zdict['xmax'] = 2.0 if zdict['x'] == 'w' else 1.0 # For a1Max = a2Max = 1.0 and R1 = R2 = 2.0.
                zdict['ymax'] = 2.0 if zdict['y'] == 'w' else 1.0 # For a1Max = a2Max = 1.0 and R1 = R2 = 2.0.
                zdict['xlabel'] = dfz.loc[dfz.z == zdict['x'], 'title'].values[0]
                zdict['ylabel'] = dfz.loc[dfz.z == zdict['y'], 'title'].values[0]

        return self

    def chart(self, dfts):

        fig, axs = plt.subplots(nrows=len(zdictss), ncols=len(zdictss[0]), figsize=(width, height), sharex=True, constrained_layout=False, squeeze=False)

        if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        for row, (rowax, zdicts) in enumerate(zip(axs, zdictss)):
            for ax, zdict in zip(rowax, zdicts):
                df = dfts[zdict['treatment']]
                x = df.loc[df[xname] == self.x_value, zdict['x']]
                y = df.loc[df[yname] == self.y_value, zdict['y']]
                ax.scatter(x, y, alpha=0.1, s=5.0)
                ax.set_xlim(0.0, zdict['xmax'])
                ax.set_ylim(0.0, zdict['ymax'])
                ax.set_xlabel(zdict['xlabel'])
                ax.set_ylabel(zdict['ylabel'])
                ax.tick_params(axis='x', labelsize=fstick)
                ax.tick_params(axis='y', labelsize=fstick)
                ax.set_box_aspect(1)

        plt.savefig(outfile, dpi=100)
        plt.close()

class Scatter:

    def prepare(self, dfs):
        for zdicts in zdictss:
            for zdict in zdicts:
                zdict['bubble_size'] = dfz.loc[dfz.z == zdict['name'], 'bubble_size'].values[0]

    def chart(self, dfts):

        fig, axs = plt.subplots(nrows=len(zdictss), ncols=len(zdictss[0]), figsize=(width, height), sharex=True, sharey=True, constrained_layout=False, squeeze=False)
        fig.supxlabel(t=x_label, y=0.02, fontsize=fslabel)
        fig.supylabel(t=y_label, x=0.04, fontsize=fslabel, ha='center')

        if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

        [dfts[d].sort_values(by=[xname, yname], inplace=True) for d in module.dirs]

        for row, (rowax, zdicts) in enumerate(zip(axs, zdictss)):
            for ax, zdict in zip(rowax, zdicts):
                if row == 0: ax.set_title(zdict['title'], pad=10.0, fontsize=fstitle)
                ax.tick_params(axis='x', labelsize=fstick)
                ax.tick_params(axis='y', labelsize=fstick)
                if x_log: ax.set_xscale('log', base=2)
                if y_log: ax.set_yscale('log', base=2)
                df = dfts[zdict['treatment']]
                x = df[xname]
                y = df[yname]
                s = df[zdict['name']]
                difs = dfts[zdict['control']][zdict['name']]/s
                if ('Grain' in zdict['name']) or ('BD' in zdict['name']): difs = 1.0/difs
                color = []
                [color.append(dif_color(dif)) for dif in difs]
                for suffix, suffixalpha in zip(['SD', ''], [0.2, 1.0]):
                    size = s + df[zdict['name'] + suffix] if suffix == 'SD' else s
                    ax.scatter(x, y, c=color, ec=color, alpha=suffixalpha, s=size*zdict['bubble_size'])
                    ax.set_xlim(x_min, x_max)
                    ax.set_ylim(y_min, y_max)
                    ax.set_box_aspect(1)

        plt.savefig(outfile, transparent=False)
        plt.close()

def dif_color(dif):
    color = (red*pow(dif,1.7), green*pow(dif,0.7), blue*pow(dif,1.7)) if dif <= 1.0 else (red*pow(dif,-1.7), green*pow(dif,-1.7), blue*pow(dif,-0.7))
    return color

def create_figure(t):
    dfts = {}
    for d in module.dirs:
        dfts[d] = dfs[d].loc[dfs[d]['Time'] == t].copy()
    pr.chart(dfts)

zdictss = module.z

dfs = {}
for d in module.dirs:
    dfs[d] = pd.concat(map(pd.read_csv, glob.glob(os.path.join(d, '*.ics'))), ignore_index=True)

if module.ftype == 'barsone': pr = BarsOne()
elif module.ftype == 'barsall': pr = BarsAll()
elif module.ftype == 'scatter': pr = Scatter()
else:
    print('No such ftype')
    exit(1)

pr.prepare(dfs)

#laststep = dfs[module.dirs[0]].Time.iat[-1]

if module.movie:
    frames = []
    for t in dfs[module.dirs[0]].Time.unique():
        outfile = f'delete{t}.png'
        create_figure(t)
        frames.append(iio.imread(outfile))
        os.remove(outfile)
        percent = t*100/laststep
        print(f'Created {percent:.1f}% of frames', end='\r')
    print('\nAdding frames to movie...')
    iio.mimsave(module.filename + '.gif', frames)
else:
    outfile = module.filename + '.png'
    pr.chart(dfs)
#    create_figure()

end_time = time.perf_counter ()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
