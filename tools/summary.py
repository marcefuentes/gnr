#! /usr/bin/env python

import glob
import pandas as pd

df = pd.concat(map(pd.read_csv, glob.glob('*.csv')), ignore_index=True)
t = df.Time.iat[-1]
dft = df.loc[df['Time'] == t]
sort = dft.sort_values([dft.columns[0], dft.columns[1], dft.columns[2]])
sort.to_csv('summary.txt', index=False)
