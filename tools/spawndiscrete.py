#! /usr/bin/env python

from glob import glob
import os
import numpy as np
import pandas as pd

folders = ['none', 'given0']
given = 0.95
num = 21
alphas = np.linspace(0.1, 0.9, num)
logess = np.linspace(-5.0, 5.0, num)
c = 101

dfs = np.empty(2, dtype=object)
for i, folder in enumerate(folders):
    filelist = glob(os.path.join(folder, '*.csv'))
    dfs[i] = pd.concat(map(pd.read_csv, filelist),
                       ignore_index=True)

df = dfs[0]
ts = df.Time.unique()
t = ts[-1]
Zs = np.empty((2, 21, 21), dtype=np.float32)

for i, folder in enumerate(folders):
    df = dfs[i]
    m = df.Time == t
    df = df.loc[m]
    Z = pd.pivot_table(df,
                       values='a2Seenmean',
                       index='alpha',
                       columns=['logES'])
    Z = Z.sort_index(axis=0, ascending=True)
    Zs[i] = Z.to_numpy()

for a, alpha in enumerate(alphas):
    for l, loges in enumerate(logess):
        filename = str(c) + '.glo'
        f = open(filename, 'w')

        f.write('Seed,1\n')
        f.write('N,12\n')
        f.write('Runs,30\n')
        f.write('Time,20\n')
        f.write('Periods,3\n')
        f.write('a1Max,1.0\n')
        f.write('a2Max,1.0\n')
        f.write('R1,2\n')
        f.write('R2,2\n')
        f.write('a2Init,0.5\n')
        f.write('ChooseGrainInit,1.0\n')
        f.write('MimicGrainInit,1.0\n')
        f.write('a2MutationSize,-6\n')
        f.write('GrainMutationSize,-6\n')
        f.write('DeathRate,-7\n')
        f.write('GroupSize,2\n')
        f.write('ChooseCost,-14\n')
        f.write('MimicCost,-14\n')
        f.write('PartnerChoice,0\n')
        f.write('Reciprocity,0\n')
        f.write('Discrete,1\n')
        f.write(f'a2low,{Zs[0, a, l]}\n')
        f.write(f'a2high,{Zs[1, a, l]}\n')
        f.write('IndirectR,0\n')
        f.write(f'alpha,{alpha:.6}\n')
        f.write(f'logES,{loges}\n')
        f.write(f'Given,{given}\n')

        f.close()
        c = c + 1
