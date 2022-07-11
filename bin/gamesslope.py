#!/usr/bin/env python

import math
import numpy as np
import matplotlib as plt
import pandas as pd

# Parameters

alpha = 0.5
R1 = 2.0
R2 = 2.0
deathrate = 1.0



width = 10.0
height = 6.0

red = 0.97
green = 0.97
blue = 0.97

fs = 14

ess = np.linspace(-5, 5, num=11)
ess = pow(2, ess)
rhos = 1.0 - 1.0/ess
givens = np.linspace(1.0, 0.0, num=11)
aC = 0.5
aD = 0.2
q1C = R1*(1.0-aC);
q1D = R1*(1.0-aD);

def fitness(q1, q2, rho):
    if rho == 0.0:
        w = pow(q1, alpha)*pow(q2, 1.0 - alpha)
    else:
        w = pow(alpha*pow(q1, rho) + (1.0 - alpha)*pow(q2, rho), 1.0/rho)
    return w

diff0 = []
diff1 = []

for rho in rhos:
    R = fitness(q1C, R2*aC, rho)
    P = fitness(q1D, R2*aD, rho)
    for given in givens:
        T = fitness(q1D, R2*(aD*(1.0-given) + aC*given), rho)
        S = fitness(q1C, R2*(aC*(1.0-given) + aD*given), rho)
        wC0 = S
        wD0 = P
        wC1 = R
        wD1 = T

        diff0.append([1/(1-rho), given, 1, abs((P-S)*0.00001/(R+P-S-T+0.0000001)), 0.0])
        diff1.append([1/(1-rho), given, 1, 0.0, 0.0])

df = pd.DataFrame(diff0, columns=['ES', 'Given', 'Time', 'wmedian', 'wmedianSD'])
df = pd.DataFrame(diff1, columns=['ES', 'Given', 'Time', 'wmedian', 'wmedianSD'])


fig = plt.subplots(nrows=len(module.rows), ncols=len(module.traits), figsize=(width, height+1.0), sharex=True, sharey=True, constrained_layout=False, squeeze=False)
fig.supxlabel(t=dfglos.loc[module.glos['x'], 'label'], y=0.02, fontsize=fslabel)
fig.supylabel(t=dfglos.loc[module.glos['y'], 'label'], x=0.04, fontsize=fslabel, ha='center')

if module.movie: fig.text(0.93, 0.02, f'Time = {t}', fontsize=14, color='grey', ha='right')

[dfts[folder].sort_values(by=[module.glos['x'], module.glos['y']], inplace=True) for folder in folderlist]

for rowax, row in zip(axs, module.rows):
    for ax, trait, traitfolder in zip(rowax, module.traits, row):
        df = dfts[traitfolder['treatment']]
        x = df[module.glos['x']]
        y = df[module.glos['y']]
        if row == module.top_row: ax.set_title(dftraits.loc[trait, 'label'], pad=10.0, fontsize=fstitle)
        ax.tick_params(axis='x', labelsize=fstick)
        ax.tick_params(axis='y', labelsize=fstick)
        if dfglos.loc[module.glos['x'], 'log']: ax.set_xscale('log', base=2)
        if dfglos.loc[module.glos['y'], 'log']: ax.set_yscale('log', base=2)
        s = df[trait]
        difs = dfts[traitfolder['control']][trait]/s
        if ('Grain' in trait) or ('BD' in trait): difs = 1.0/difs
        if '6' in trait: difs = 1.0/difs
        color = []
        [color.append(dif_color(dif)) for dif in difs]
        for suffix, suffixalpha in zip(['SD', ''], [0.2, 1.0]):
            size = s + df[trait + suffix] if suffix == 'SD' else s
            ax.scatter(x, y, c=color, ec=color, alpha=suffixalpha, s=size*dftraits.loc[trait, 'bubble_size'])
            ax.set_xlim(self.glosx_min, self.glosx_max)
            ax.set_ylim(self.glosy_min, self.glosy_max)
            ax.set_box_aspect(1)

plt.savefig(outfile, transparent=False)
plt.close()

