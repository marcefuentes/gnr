#! /usr/bin/env python

from glob import glob
import os
import time

import numpy as np
import pandas as pd

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

folders = ['given100', 'given95', 'given50']
subfolders = ['none', 'p', 'r']

# Get data

def read_file(file, alltimes):
    df = pd.read_csv(file)
    if not alltimes:
        df = df.tail(1)
    return df

filelist = glob('given00/none/*.csv')
df = [read_file(file, False) for file in filelist]
df = pd.concat(df, ignore_index=True)
df = df.sort_values(['alpha', 'logES'])
high = df['a2Seenmean'].values
high = np.tile(high, (len(folders), 1)).flatten()
print(high.shape)

dfs = []
for folder in folders:
    df_list = []
    for subfolder in subfolders:
        filelist = glob(os.path.join(folder, subfolder, '*.csv'))
        df = [read_file(file, False) for file in filelist]
        df_concat = pd.concat(df, ignore_index=True)
        df_list.append(df_concat)
    dfs.append(df_list)

dfns = []
for i, subfolder in enumerate(subfolders):
    df_concat = pd.concat([dfs[j][i] for j in range(len(folders))])
    dfns.append(df_concat)

for df in dfns:
    df = df.sort_values(['Given', 'alpha', 'logES'])

df = dfns[0]
low = df['a2Seenmean'].values
given = df['Given'].values
alpha = df['alpha'].values
loges = df['logES'].values
rho = 1.0 - 1.0/pow(2.0, loges)
print(low.shape)

T = my.fitness(high, low, given, alpha, rho)
R = my.fitness(high, high, given, alpha, rho)
P = my.fitness(low, low, given, alpha, rho)
S = my.fitness(low, high, given, alpha, rho)

X = np.stack((T, R, P, S), axis=1)
print(X.shape)

CG = dfns[1]['ChooseGrainmean'].values
MG = dfns[2]['MimicGrainmean'].values
print(CG.shape)
print(MG.shape)

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
