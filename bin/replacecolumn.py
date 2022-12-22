#!/usr/bin/python

import os
import pandas as pd

alpha = 0.75

for root, dirs, files in os.walk('.'):
    for filename in files:
        if ('csv' in filename) or ('frq' in filename):
            df = pd.read_csv(os.path.join(root, filename))
            newcol = [alpha]*len(df.index)
            df['alpha'] = newcol
            df.drop('N', axis=1, inplace=True)
            df.insert(2, 'alpha', df.pop('alpha'))
            df.to_csv(os.path.join(root, filename), float_format='%.6f', index=False)
